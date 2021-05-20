#Calling libraries
import argparse
import cosima_cookbook as cc
import netCDF4 as nc
import xarray as xr
import numpy as np
import pandas as pd
import copy
import os
import re
import rasterio
import geopandas
import rasterio.plot
import rioxarray
from shapely.geometry import mapping, Polygon
import calendar

########
#Defining functions

########
#Loads ACCESS-OM2-01 data for the Southern Ocean
def getACCESSdata(var, start, end, freq, ses, minlat = -90, maxlat = -45, exp = '01deg_jra55v140_iaf_cycle2'):
    '''
    Inputs:
    var - Short name for the variable of interest
    start - Time from when data has to be returned
    end - Time until when data has to be returned
    freq - Time frequency of the data
    ses - Cookbook session
    minlat - minimum latitude from which to return data. If not set, defaults to -90 to cover the Southern Ocean.
    maxlat - maximum latitude from which to return data. If not set, defaults to -45 to cover the Southern Ocean.
    exp - Experiment name
        
    Output:
    Data array with corrected time and coordinates within the specified time period and spatial bounding box.
    '''
    #Accessing data
    vararray = cc.querying.getvar(exp, var, ses, frequency = freq, start_time = start, end_time = end)
    #Subsetting data to area of interest
    vararray = vararray.sel(yt_ocean = slice(minlat, maxlat))
    return vararray

########
#Correcting longitude values in a data array so they are between -180 and +180 degrees
def corrlong(array):
    '''
    Inputs:
    array - Data array on which longitude corrections will be applied.
    
    Output:
    Data array with corrected longitude values.
    '''
    
    #Making a deep copy of original longitude values in the array being corrected
    corr_lon = copy.deepcopy(array.xt_ocean.values)
    
    ##Now we need to correct any values smaller than -180 and replace them with values between +80 and +180. Note that the smallest longitude value (-279.95) should be +80.05.
    #While -180.05 should have a correct value of +179.95.
    corr_lon[np.where(corr_lon < -180)] = sorted(-corr_lon[np.where((corr_lon >= -180) & (corr_lon <= -80))])
    
    #Assign corrected longitude values to the data array being corrected
    array.coords['xt_ocean'] = corr_lon
    
    #Longitude values must be sorted from smallest to largest prior to plotting
    array = array.sortby(array.xt_ocean)
    
    return array

########
#This function assigns the same coordinate reference system (CRS) to the data array being clipped as the clipping shapefile and then clips the data array
def clipDataArray(array, shp):
    '''
    Inputs:
    array - Data array to be clipped.
    shp - Shapefile to be used for clipping.
    
    Output:
    Clipped data array.
    '''
    crs = 'epsg:' + str(shp.crs.to_epsg())
    array.rio.set_spatial_dims(x_dim = 'xt_ocean', y_dim = 'yt_ocean', inplace = True)
    array.rio.write_crs(crs, inplace = True)
    
    clipped = array.rio.clip(shp.geometry.apply(mapping), shp.crs, drop = True, invert = False, all_touched = False)
    
    return clipped

########
#Calculates weighted means by season, by month or per timestep
def weightedMeans(array, weights, meanby = 'timestep'):
    '''
    Inputs:
    array - Data array containing variable from which means will be calculated
    weights - Data array containing weights
    meanby - Define how means will be calculate: timestep, month or season. Default set to 'timestep'
    
    Output:
    Data array containing weighted means
    '''
      
    #Calculate weights
    weights = weights/weights.sum()
        
    #Apply weights to variable - Calculate weighted mean over timestep and then calculate the mean by season
    if meanby == 'season':
        weighted_mean = (array*weights).groupby('time').sum(('xt_ocean', 'yt_ocean')).groupby('time.season').mean()
    elif meanby == 'month':
        weighted_mean = (array*weights).groupby('time').sum(('xt_ocean', 'yt_ocean')).groupby('time.month').mean()
    elif meanby == 'timestep':
        weighted_mean = (array*weights).groupby('time').sum(('xt_ocean', 'yt_ocean'))
        
    return weighted_mean

########
#This function adds a time dimension containing the year of sampling in an array containing summarised data. It must be specified is the data is summarised by season or per month
def addTimeYear(array, year, by = 'season'):
    '''
    Inputs:
    array - Data array to which a time dimension containing the year of sampling will be added to an array that contains four timesteps, one per season
    year - A string containing the year that will be added as a time dimension
    by - A string describing if data array is group by months (monthly), or per season (season). Default set to season
    
    Output:
    Data array with a time dimension containing the year of sampling
    '''
    if not isinstance(year, str):
        year = str(year)
    
    if by == 'season':
        #Create a time variable to add as a coordinate to each array
        time = [year]*4
        #Add time coordinate to data array
        x = array.assign_coords(time = ('season', time))
    elif by == 'month':
        time = [year]*12
        #Add time coordinate to data array
        x = array.assign_coords(time = ('month', time))
    
    #Return data array with the time dimension
    return x

########
#This function access all netcdf files that are contained in one folder and have a particular keyword in their name
def stackData(folder, keyword):
    '''
    Inputs:
    folder - The location of the folder where all netcdf files to be concatenated are located
    keyword - String containing a keyword that will be used to identify all files to be concatenated.
    
    Output:
    Concatenated data array
    '''
    
    if not isinstance(keyword, str):
        print('Keyword argument must be a string.')
    
    #Get a list files that contain the keyword provided and order them alphabetically
    filelist = sorted(list(filter(re.compile('.*' + keyword + '.*').match, os.listdir(folder))))
    
    #Create an empty list that will contain calculations for every year
    combData = []

    #Loop for every year included in the analysis
    for i in np.arange(0, len(filelist)):
        #Read file
        x = xr.open_dataarray(os.path.join(folder, filelist[i]))
        combData.append(x)

    #Concatenate all items included in the list created in the loop to create one data array per sector
    comb_data = xr.concat(combData, dim = 'season')

    #Return concatenated data array
    return comb_data


########
#This function corrects year values
def corrYears(xarray):
    Months = [calendar.month_abbr[m] for m in xarray.month.values]
    xarray.coords['season'] = xarray['time'].values[:,0]
    xarray.coords['month'] = Months
    
    
########
#This function can be used to combine various data arrays into one
def getFileList(filein, yrs):
    '''
    Inputs:
    filein - refers to the file path of the folder containing the datasets to be used in calculations.
    yrs - is a numpy array containing a list of years of interest for calculations.
      
    Outputs:
    Three lists: adv_list, ret_list, and sea_list which contain the list of files containing sea ice advance, retreat and total season duration respectively.
    '''
    #List netcdf files containing sea ice seasonality data
    filelist = os.listdir(filein)
 
    #Extract files for baseline years
    basefiles = []
    for i in filelist:
        for j in yrs:
            if i[10:14] == str(j):
                basefiles.append(i)
    #Remove variables no longer in use
    del filelist

    #Separate files based on whether they contain information about sea ice advance, retreat or season. Order them alphabetically.
    adv_list = sorted([i for i in basefiles if 'Adv' in i])
    ret_list = sorted([i for i in basefiles if 'Ret' in i])
    sea_list = sorted([i for i in basefiles if 'Dur' in i])
    #Remove variables no longer in use
    del basefiles
    
    #Return file lists
    return adv_list, ret_list, sea_list


########
#This function can be used to calculate baseline means or to extract files for any other time period
def combineData(filepath, filelist, dir_out):
    '''
    Inputs:
    filepath - refers to the file path of the folder containing the netcdf files to be combined into a single data array.
    filelist - contains the actual names of the netcdf files that will be combined..
    dir_out - file path of the folder where combined data arrays will be saved
      
    Outputs:
    Three dimensional data array containing all files provided in the filelist input. The data array is saved to the path provided in dir_out and it can also be assigned to a variable.
    '''
    #Create variable to hold combined data arrays
    combData = []
    #Create loop based on the total length of filelist
    for i in np.arange(0, len(filelist)):
        #Open data array
        x = xr.open_dataarray(os.path.join(filepath, filelist[i]), autoclose = True)
        #Put all files in one variable
        combData.append(x)
        
    #Create one data array with the data contained in the combined variable
    combined = xr.concat(combData, dim = 'time')
    combined.to_netcdf(os.path.join(dir_out, (filelist[0][0:15]+filelist[-1][15:19]+'.nc')))
    return combined
    
########
def main(inargs):
    '''Run the program.'''

if __name__ == '__main__':
    description = 'This script contains functions used to perform timeseries within different sectors of the Southern Ocean.'
    parser = argparse.ArgumentParser(description = description)

    args = parser.parse_args()
    main(args)