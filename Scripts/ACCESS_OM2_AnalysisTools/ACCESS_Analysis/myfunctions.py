#Calling libraries
import cosima_cookbook as cc
import xarray as xr


#Defining functions
########
def getACCESSdata(var, start, end, freq, ses, minlat = -90, maxlat = -45, exp = '01deg_jra55v140_iaf_cycle2', ice_data = False):
    '''
    Loads ACCESS-OM2-01 sea ice and ocean data for the Southern Ocean. If ice data is accessed, it corrects the time and coordinate grid to match ocean outputs.
    
    :param var: Variable of interest
    :type var: str
    
    :param start: Date in the form YYYY, YYYY-MM or YYYY-MM-DD from when the data must be returned
    :type start: str
    
    :param end: Date in the form YYYY, YYYY-MM or YYYY-MM-DD until when the data must be returned
    :type end: str
    
    :param freq: Frequency of the data, daily or monthly
    :type freq: str
    
    :param ses: COSIMA cookbook session
    :type ses: sqlalchemy.orm.session.Session
    
    :par minlat: minimum latitude from which to return data. If not set, defaults to -90 to cover the Southern Ocean.
    :type minlat: int
    
    :par maxlat: maximum latitude from which to return data. If not set, defaults to -45 to cover the Southern Ocean.
    :type maxlat: int
    
    :par exp: Experiment name. Default is 01deg_jra55v140_iaf_cycle2.
    :type exp: str
    
    :par ice_data: When True the variable being called is related to sea ice, when False is not. Default is set to False (i.e., it assumes variable is related to the ocean).
    :type ice_data: bool
        
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