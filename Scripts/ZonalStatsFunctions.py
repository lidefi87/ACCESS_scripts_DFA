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
import statsmodels.api as sm
import datetime as dt
import scipy.stats as ss

########
#Defining functions

########
#Loads ACCESS-OM2-01 sea ice and ocean data for the Southern Ocean. If ice data is accessed, it corrects the time and coordinate grid to match ocean outputs.
def getACCESSdata(var, start, end, freq, ses, minlat = -90, maxlat = -45, 
                  exp = '01deg_jra55v140_iaf_cycle2', ice_data = False):
    '''
    Inputs:
    var - Short name for the variable of interest
    start - Time from when data has to be returned
    end - Time until when data has to be returned
    freq - Time frequency of the data
    ses - Cookbook session
    minlat - minimum latitude from which to return data. If not set, defaults to -90 to cover the Southern Ocean.
    maxlat - maximum latitude from which to return data. If not set, defaults to -45 to cover the Southern Ocean.
    exp - Experiment name. Default is 01deg_jra55v140_iaf_cycle2.
    ice_data - Boolean, when True the variable being called is related to sea ice, when False is not. Default is set to False (i.e., it assumes variable is related to the ocean).
        
    Output:
    Data array with corrected time and coordinates within the specified time period and spatial bounding box.
    '''
    #Accessing data
    vararray = cc.querying.getvar(exp, var, ses, frequency = freq, start_time = start, end_time = end)
    
    #If data being accessed is an ice related variable, then apply the following steps
    if ice_data == True:
        #Accessing corrected coordinate data to update geographical coordinates in the array of interest
        geolon_t = cc.querying.getvar(exp, 'geolon_t', ses, n = -1)
        #Apply time correction so data appears in the middle (12:00) of the day rather than at the beginning of the day (00:00)
        vararray['time'] = vararray.time.to_pandas() - dt.timedelta(hours = 12)
        #Change coordinates so they match ocean dimensions 
        vararray.coords['ni'] = geolon_t['xt_ocean'].values
        vararray.coords['nj'] = geolon_t['yt_ocean'].values
        #Rename coordinate variables so they match ocean data
        vararray = vararray.rename(({'ni':'xt_ocean', 'nj':'yt_ocean'}))
    
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
        
    #Set the spatial dimensions of the xarray being clipped
    array.rio.set_spatial_dims(x_dim = 'xt_ocean', y_dim = 'yt_ocean', inplace = True) #inplace = True updates the array instead of creating a copy
    #Assign a CRS to the xarray that matches the shapefile used for clipping. CRS included is CF compliant.
    array.rio.write_crs(shp.crs, inplace = True) #inplace = True updates the array instead of creating a copy
    
    #Clipping maintains only those pixels whose center is within the polygon boundaries and drops any data not meeting this requirement.
    clipped = array.rio.clip(shp.geometry, shp.crs, drop = True, invert = False, all_touched = False)
    
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
    filelist - contains the actual names of the netcdf files that will be combined.
    dir_out - file path of the folder where combined data arrays will be saved.
      
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
#This function creates a colour palette using Crameri's palettes (Crameri, F. (2018), Scientific colour-maps, Zenodo, doi:10.5281/zenodo.1243862)
def colourMaps(colourLibraryPath, palette, rev = True):
    '''
    Inputs:
    colourLibraryPath - the file path where the palettes are currently saved.
    palette - name of the palette to be created.
    rev - Boolean. If True, it will create a normal and reversed version of the palette. If False, it will only return one palette
    
    Outputs:
    One or two palettes based on Crameri (2018) that can be used to colour matplotlib figures
    '''
    #Load relevant libraries to set Scientific Colour Map library
    from matplotlib.colors import LinearSegmentedColormap
    from matplotlib.colors import ListedColormap

    #Set path where the scientific library is found
    cm_data = np.loadtxt(os.path.join(colourLibraryPath, palette, (palette + '.txt')))
    #Create a colour map based on 'palette' argument
    pal_map_adv = LinearSegmentedColormap.from_list(palette, cm_data)
        
    if rev == True:
        pal_map_ret = ListedColormap(cm_data[::-1])
        return pal_map_adv,pal_map_ret
    else:
        return pal_map_adv
    
    
########
#This function performs a linear trend calculation and returns the coefficients as well as p-values for the linear regression
def linearTrends(y, x, rsquared = False):
    '''
    Inputs:
    y - data array with information about dependent variable
    x - data array with information about independent variable
    rsquared - Boolean. If set to True then r squared values will be calculated and returned as outputs
        
    Output:
    Coefficients and p-values of linear regression
    '''
    #To check extra information available in the model use
        #dir(model.fit())
        
    if rsquared == True:
        model = sm.OLS(y, x)
        coef = model.fit().params[1]
        sig = model.fit().pvalues[1]
        rsq_adj = model.fit().rsquared_adj
        return coef, sig, rsq_adj
    else:
        model = sm.OLS(y, x)
        coef = model.fit().params[1]
        sig = model.fit().pvalues[1]
        return coef, sig

    
########
#Defining function that will be applied across dimensions
def lm_yr(y, x):
    '''
    Inputs:
    x - refers to the independent variable
    y - refers to the dependent variable
    dim - refers to the dimension along which the function will be mapped
    
    Output:
    Data array containing the results of a simple linear model
    '''
    return ss.linregress(x, y)
    
    
########
#This function calculates anomalies 
def AnomCalc(array, clim_array, std_anom = False):
    '''
    Inputs:
    array - refers to a data array containing information for the period being compared to the baseline. It could include just one year or multiple years (decades)
    clim_array - three dimensional array containing data over the baseline period from which anomalies will be calculated
    std_anom - boolean variable that if set to True will result in standarised anomalies being calculated
      
    Outputs:
    Data array containing anomalies.
    '''
    
    #Calculate long term mean of array
    m_climarray = clim_array.mean('time')
      
    #Calculate anomalies
    #Standarised anomalies
    if std_anom == True:
        s_climarray = clim_array.std('time')
        anom = (array - m_climarray)/s_climarray
    #Non-standarised anomalies
    elif std_anom == False:
        anom = array - m_climarray
    
    #Return anomalies
    return anom


########
#Getting maximum and minimum years included in the analysis
def colbarRange(dict_data, sector, season):
    '''
    Inputs:
    dict_data - refers to a dictionary that contains the datasets being plotted.
    sector - refers to a list of sectors
    season - refers to a list of seasons
    '''
    
    #Define variables that will store the maximum and minimum values
    maxV = []
    minV = []
    
    sector = sector*len(season)
    season = np.concatenate([[i]*len(sector) for i in season], axis = 0)
    for sec, sea in zip(sector, season):
        maxV.append(dict_data[f'{sec}_{sea}'].indexes['time'].year.max())
        minV.append(dict_data[f'{sec}_{sea}'].indexes['time'].year.min())
        
    maxV = max(maxV)
    minV = min(minV)

    return minV, maxV


########
def SeaIceAdvArrays(array, thres = 0.15, ndays = 5, **kwargs):
    '''
    The `SeaIceAdvArrays` function below was losely based on the `calc_ice_season` function from the `aceecostats` R package developed by Michael Sumner at AAD. This function calculates annual sea ice advance, retreat and total sea ice season duration as defined by Massom et al 2013 [DOI:10.1371/journal.pone.0064756].  
Briefly, if sea ice concentration in any pixel is at least 15% over five consecutive days, sea ice is considered to be advancing. Sea ice is retreating when its concentration is below 15% in any pixel until the end of the sea ice year. Sea ice season duration is the period between day of advance and retreat. Sea ice year is between February 15 and February 14 the following year.
    Inputs:
    array is the data array on which sea ice seasonality calculations will be performed
    dir_out is the file path to the folder where outputs should be saved.
    thres refers to the minimum sea ice concentration threshold. The default is set to 0.15
    ndays is the minimum amount of consecutive days sea ice must be above threshold to be classified as advancing. Default set to 5
    
    Outputs:
    Function saves three data arrays as netcdf files: advance, retreat and season duration. Data arrays can also be saved as variables in the notebook.
    '''
    
    #Extracting maximum and minimum year information to extract data for the sea ice year 
    MinY = str(array.time.dt.year.values.min())
    MaxY = str(array.time.dt.year.values.max())
    
    #Selecting data between Feb 15 to Feb 14 (sea ice year)
    array = array.sel(time = slice(f'{MinY}-02-15', f'{MaxY}-02-14'))
    
    ########
    #Preparing masks to perform calculations on areas of interest only
    #Calculate timesteps in dataset (365 or 366 depending on whether it was a leap year or not)
    timesteps = len(array.time.values)

    #Identify pixels (or cells) where sea ice concentration values are equal or above the threshold
    #Resulting data array is boolean. If condition is met then pixel is set to True, otherwise set to False
    threshold = xr.where(array >= thres, True, False)

    #Creating masks based on time over threshold
    #Add values through time to get total of days with ice cover of at least 15% within a pixel
    rsum = threshold.sum('time')

    #Boolean data arrays for masking
    #If the total sum is zero, then set pixel to be True, otherwise set to False. 
    #This identifies pixels where minimum sea ice concentration was never reached.
    noIce = xr.where(rsum == 0, True, False)
    #If the total sum is less than the minimum days, then set pixel to be True, otherwise set to False. 
    #This identifies pixels where sea ice coverage did not meet the minimum consecutive days requirement.
    noIceAdv = xr.where(rsum < ndays, True, False)
    #If the total sum is the same as the timesteps, then set pixel to be True, otherwise set to False.
    #This identifies pixels where sea ice concentration was always at least 15%
    alwaysIce = xr.where(rsum == timesteps, True, False)
    #Remove unused variables
    del rsum

    ########
    #Sea ice advance calculations
    #Use cumulative sums based on time. If pixel has sea ice cover below threshold, 
    #then cumulative sum is reset to zero
    adv = threshold.cumsum(dim = 'time')-threshold.cumsum(dim = 'time').\
    where(threshold == 0).ffill(dim = 'time').fillna(0)
    #Note: ffill adds nan values forward over a specific dimension

    #Find timestep (date) where the minimum consecutive sea ice concentration was first 
    #detected for each pixel
    #Change all pixels that do not meet the minimum consecutive sea ice concentration to False. 
    #Otherwise maintain their value.
    advDate = xr.where(adv == ndays, adv, False)
    #Find the time step index where condition above was met.
    advDate = advDate.argmax(dim = 'time')
    #Apply masks of no sea ice advance and sea ice always present.
    advDate = advDate.where(noIceAdv == False, np.nan).where(alwaysIce == False, 1)
    #Remove unused variables
    del adv

    ########
    #Sea ice retreat calculations
    #Reverse threshold data array (time wise) - So end date is now the start date and calculate 
    #cumulative sum over time
    ret = threshold[::-1].cumsum('time')
    del threshold
    #Change zero values to 9999 so they are ignored in the next step of our calculation
    ret = xr.where(ret == 0, 9999, ret)
    #Find the time step index where sea ice concentration change to above threshold.
    retDate = ret.argmin(dim = 'time')
    #Substract index from total time length
    retDate = timesteps-retDate
    #Apply masks of no sea ice over threshold and sea ice always over threshold.
    retDate = retDate.where(noIce == False, np.nan).where(alwaysIce == False, timesteps)
    #Remove unused variables
    del ret
    
    ########
    #Sea ice duration
    durDays = retDate-advDate
    #Remove unused variables
    del noIce, noIceAdv, alwaysIce
    
    ########
    #Adding a time dimension to newly created arrays and removing unused dimensions
    def addTime(array, year):
        #Create a time variable to add as dimension to each array - Only one timestep included
        time = pd.date_range(f'{year}-02-15', f'{year}-02-16', freq = 'D', closed = 'left')
        #Add time dimension to data array
        x = array.expand_dims({'time': time}).assign_coords({'time': time})
        #Remove dimensions that are not needed
        x = x.drop('TLON').drop('TLAT').drop('ULON').drop('ULAT')
        #Return 
        return x
         
    #Applying function
    advDate2 = addTime(advDate, MinY)
    retDate2 = addTime(retDate, MinY)
    durDate = addTime(durDays, MinY)
    del advDate, retDate, durDays
      
    ########
    #Save corrected outputs as netcdfiles
    if 'dir_out' in kwargs.keys():
        #Check output folder exists
        os.makedirs(kwargs.get('dir_out'), exist_ok = True)
    
        #Define output paths
        advpath = os.path.join(dir_out, (f'SeaIceAdv_{MinY}-{MaxY}.nc'))
        retpath = os.path.join(dir_out, (f'SeaIceRet_{MinY}-{MaxY}.nc'))
        durpath = os.path.join(dir_out, (f'SeaIceDur_{MinY}-{MaxY}.nc'))
    
        #Save files simultaneously
        xr.save_mfdataset(datasets = [advDate2.to_dataset(), 
                                      retDate2.to_dataset(), 
                                      durDate.to_dataset()], 
                          paths = [advpath, retpath, durpath])
    
    #Return data arrays as outputs
    return (advDate2, retDate2, durDate)


########
def main(inargs):
    '''Run the program.'''

if __name__ == '__main__':
    description = 'This script contains functions used to perform timeseries within different sectors of the Southern Ocean.'
    parser = argparse.ArgumentParser(description = description)

    args = parser.parse_args()
    main(args)