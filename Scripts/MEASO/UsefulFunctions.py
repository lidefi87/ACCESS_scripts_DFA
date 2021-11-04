#Calling libraries
import argparse
import netCDF4 as nc
import xarray as xr
import numpy as np
import pandas as pd
import copy
import os
import re
import calendar
import statsmodels.api as sm
from clef.code import *
from glob import glob

########
#Defining functions

########
#Search data within clef database
def searchACCESS(var, model, freq, exp, **kwargs):
    '''
    Inputs:
    var - str, code for variable of interest used in CMIP6 models
    model - str, CMIP6 model where variable of interest will be search
    freq - str, requency of the data needed
    
    Optional inputs:
    variant - str, experiment run to be search
    time_frame - list, time period of interest given to the nearest decade 
    (e.g., if interested in the period between 1991 and 2008, then use: [1990, 2000, 2010]
    
    Returns:
    files - list, includes file paths for all datasets that meet search requirements
    within the specified time period
    '''
    
    #Creating a session and connecting to database
    db = connect()
    s = Session()
    
    #Extract variables from kwargs
    if 'variant' in kwargs.keys():
        variant = kwargs.get('variant')
    else:
        variant = None
    
    #Creating dictionary to perform search using clef
    search_dict = {'variable_id': var, 'model': model, 'frequency': freq, 
              'experiment_id':exp, 'variant_label': variant}
    df = search(s, project = 'CMIP6', latest = True, **search_dict)
    
    #If experiment runs (variants) need to be search, use the code below
    # sorted([re.search('r[0-9]{,2}i[0-9]{,2}p[0-9]{,2}f[0-9]{,2}', 
    #                   f).group() for f in df['path']])
    
    if 'time_frame' in kwargs.keys():
        time_frame = kwargs.get('time_frame')
        #Creating a file path list of all available datasets meeting search
        #parameters
        folder_path = [f for f in df['path'] if variant in f]
        filenames = sorted(glob(os.path.join(folder_path[0], '*.nc')))
    
        #Subsetting list to select only files within the time period of interest
        files = []
        for yr in time_frame:
            files.append([f for f in filenames if str(yr) in f][0])
    
    else:
        files = [os.path.join(df['path'][0], list(df['filename'][0])[0])]
    
    return files

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
    lon = array['longitude'][0].values
    
    #Values changed from 0-360 to between -180 and +180.
    lon[np.where(lon > 180)] = lon[np.where(lon >180)]-360
    
    #Assigning corrected longitudes and latitude values to i and j, which are index based
    array.coords['i'] = lon
    array.coords['j'] = array['latitude'][:,0].values
    
    #Removing latitude and longitude dimensions and renaming i and j dimensions
    array = array.drop(('latitude', 'longitude'))
    array = array.rename(({'i':'longitude', 'j':'latitude'}))
    
    #Longitude values must be sorted from smallest to largest prior to plotting
    array = array.sortby(array.longitude)
    
    return array

########
#This function loads the data that was found to meet search requirements when the searchACCESS function was applied
def loadData(filelist, var_name, SO = True, weights = False, **kwargs):
    '''
    Inputs:
    filelist - list, filepaths for variables of interest
    var_name - str, code for variable of interest used in CMIP6 models
    SO - boolean, if True it will return data frame for the Southern Ocean
    weights - boolena, if True it will only return one time slice to be used as weights (e.g., area, volume, height) for calculating weighted statistics
    
    Optional inputs:
    years - list, years to be included in the data
    months - list, must be a string with the start and end month as two digits
    depth_range - list, including maximum and minimum depths to be selected
        
    Returns:
    da - data frame, includes data for variable, time period and habitat of interest
    '''
    
    #Empty variable to store data frames
    var = []
    
    if len(filelist) > 1:
        #Looping through files and stacking them
        for f in filelist:
            var.append(xr.open_dataset(f, mask_and_scale = True))
        #Concatenating files across time dimension
        var = xr.concat(var, dim = 'time')
    else:
        var = xr.open_dataset(filelist[0], mask_and_scale = True)
        if weights == True:
            var = var[var_name][0]
    
    if 'years' in kwargs.keys() and 'months' in kwargs.keys():
        years = kwargs.get('years')
        months = kwargs.get('months')
        #Selecting data for the time period of interest
        s_time = f'{str(years[0])}-{str(months[0])}'
        e_time = f'{str(years[-1])}-{str(months[-1])}'
        
        var = var[var_name].sel(time = slice(s_time, e_time))
    
    if 'depth_range' in kwargs.keys():
        depths = kwargs.get('depth_range')
        #Subsetting data based on depths of interest
        var = var.sel(lev = slice(depths[0], depths[-1]))
        
    if type(var) == xr.core.dataset.Dataset:
        var = var[var_name]
        
    #Apply latitude and longitude corrections
    var = corrlong(var)
    
    if SO == True:
        #Select the SO
        var = var.sel(latitude = slice(-90, -30))
    
    #Return variable with variable of interest for specified time frame, depth and area
    return var

########
#Creates a mask from a netcdf file that can be applied to a data array
def creatingMask(mask_file):
    '''
    Inputs:
    maskfile - str, filepath for the location of mask. Must be a netcdf file
            
    Returns:
    mask_reg - data frame, to mask data
    regionNames - list, containing names of regions within mask
    '''
    
    #Loading mask
    mask = xr.load_dataarray(mask_file)
    
    #Applying mask
    #Getting region names from mask
    regionNames = sorted(set(mask.region.values))

    #Subsetting shapefiles into regions
    #Initialise dictionary that will contain sector limits
    mask_reg = {}

    #Saving each sector as an entry in the dictionary
    for reg in regionNames:
        mask_reg[f"{reg}"] = mask.sel(region = reg)

    return mask_reg, regionNames

########
#Calculates weighted means by season, by month or per timestep
def weightedMeans(regions, var_df, mask_df, weights):
    '''
    Inputs:
    regions - list, containing names of regions within mask
    var_df - data frame, containing variables of interest
    mask_df - data frame, to mask data
    weights - data frame, containing weights to be applied to mean calculations
            
    Returns:
    mean_calcs - data frame, containing weighted monthly means per sector
    '''
    #Empty lists to save results
    mean_calcs = []
    
    #For each region in mask perform calculation
    for reg in regions:
        #Apply mask to variable 
        var_reg = var_df*mask_df[reg]
        #Apply mask to volume data and replace NAs with zeroes prior to weighting data
        weight_reg = weights*mask_df[reg]
        weight_reg = weight_reg.fillna(0)
        
        #Calculate weighted means per sector
        mean_weighted_var = var_reg.weighted(weight_reg).mean(('lev', 'latitude', 'longitude'))
        #Save in empty list before concatenation
        mean_calcs.append(mean_weighted_var)
    
    # Create one netcdf file per calculation and saving result
    mean_calcs = xr.concat(mean_calcs, dim = 'region')
    
    return mean_calcs


########
#This function calculates weighted and unweighted standard deviations
def std_dev(regions, var_df, mask_df, weights, weighted_means):
    '''
    Inputs:
    regions - list, containing names of regions within mask
    var_df - data frame, containing variables of interest
    mask_df - data frame, to mask data
    weights - data frame, containing weights to be applied to mean calculations
    weighted_means - data frame, containing weighted means
            
    Returns:
    un_std_calcs - data frame, containing unweighted monthly std dev per sector
    w_std_calcs - data frame, containing weighted monthly std dev per sector
    '''
    
    #Empty lists to save results
    un_std_calcs = []
    w_std_calcs = []
    
    #For each region in mask perform calculation
    for reg in regions:
        #Apply mask to variable 
        var_reg = var_df*mask_df[reg]
        #Apply mask to volume data and replace NAs with zeroes prior to weighting data
        weight_reg = weights*mask_df[reg]
        weight_reg = weight_reg.fillna(0)
        
        #Calculating unweighted standard deviation
        std_unweighted_var = var_reg.std(('lev', 'latitude', 'longitude'))
        un_std_calcs.append(std_unweighted_var)
        
        #Calculating weighted standard deviation
        #Calculate Bessel's correction for sample std (i.e., (n-1)/n)
        bes_cor = ((weight_reg/weight_reg).sum()-1)/(weight_reg/weight_reg).sum()
        #Divisor
        divisor =  weight_reg.sum()*bes_cor
        #Dividend
        std = []
        for t in var_reg.time:
            dif = (weight_reg*((var_reg.sel(time = t)-weighted_means.sel(region = reg, 
                                                                         time = t))**2)).sum(('lev',
                                                                                              'latitude',
                                                                                              'longitude'))
            var = dif/divisor
            std.append(np.sqrt(var.values))
        std = xr.DataArray(std, dims = ['time'], coords = {'time': var_reg.time})
        std = std.expand_dims({'region': [reg]})
        w_std_calcs.append(std)
    
    # Create one netcdf file per calculation and saving result
    un_std_calcs = xr.concat(un_std_calcs, dim = 'region')
    w_std_calcs = xr.concat(w_std_calcs, dim = 'region')
    
    return un_std_calcs, w_std_calcs


########
#This function calculates percentiles
def perc_calc(regions, var_df, mask_df, percentiles):
    '''
    Inputs:
    regions - list, containing names of regions within mask
    var_df - data frame, containing variables of interest
    mask_df - data frame, to mask data
    percentiles - list, percentiles that need to be calculated
            
    Returns:
    per_calcs - data frame, containing monthly percentiles per sector
    '''
    #Empty lists to save results
    per_calcs = []

    #For each region in mask perform calculation
    for reg in regions:
        #Apply mask to variable 
        var_reg = var_df*mask_df[reg]
        
        #Calculate percentiles per sector
        perc_var = var_reg.quantile(percentiles, ('lev', 'latitude', 'longitude'))
        perc_var = perc_var.expand_dims({'region': [reg]})
        #Save in empty list before concatenation
        per_calcs.append(perc_var)
    
    #Create one netcdf file per calculation and saving result
    per_calcs = xr.concat(per_calcs, dim = 'region')
    
    return per_calcs
    
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
def linearTrends(y, x):
    '''
    Inputs:
    y - data array with information about dependent variable
    x - data array with information about independent variable
        
    Output:
    Coefficients and p-values of linear regression
    '''
    
    model = sm.OLS(y, x)
    coef = model.fit().params[1]
    sig = model.fit().pvalues[1]
    rsq_adj = model.fit().rsquared_adj
    
    #To check extra information available in the model use
    #dir(model.fit())
    
    return coef, sig, rsq_adj

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
def main(inargs):
    '''Run the program.'''

if __name__ == '__main__':
    description = 'This script contains functions used to perform timeseries within different sectors of the Southern Ocean.'
    parser = argparse.ArgumentParser(description = description)

    args = parser.parse_args()
    main(args)