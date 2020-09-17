import xarray as xr
import numpy as np

def write_stat(nc, group_name, binning):
    
    grp = nc.createGroup(group_name)
    grp.createDimension('lon', len(binning.x))
    grp.createDimension('lat', len(binning.y))
    
    longitude = grp.createVariable('lon', 'f4', 'lon', zlib=True)
    longitude[:] = binning.x
    latitude = grp.createVariable('lat', 'f4', 'lat', zlib=True)
    latitude[:] = binning.y
    
    stats = ['min', 'max', 'median', 'sum', 'sum_of_weights', 'variance', 'mean', 'count', 'kurtosis', 'skewness']
    for variable in stats:
        
        var = grp.createVariable(variable, binning.variable(variable).dtype, ('lat','lon'), zlib=True)
        var[:, :] = binning.variable(variable).T 

        
def write_timeserie_stat(data_vector, time_vector, freq, output_filename):
    
    # convert data vector and time vector into xarray.Dataarray
    da = xr.DataArray(data_vector, coords=[time_vector], dims="time")
    
    # resample 
    da_resample = da.resample(time=freq)
    
    # compute stats
    vmean = da_resample.mean()
    vminimum = da_resample.min()
    vmaximum = da_resample.max()
    vcount = da_resample.count()
    vvariance = da_resample.var()
    vmedian = da_resample.median()
    vrms = np.sqrt(np.square(da).resample(time=freq).mean())
    
    # save stat to dataset
    ds = xr.Dataset(
        {
            "mean": (("time"), vmean.values),
            "min": (("time"), vminimum.values),
            "max": (("time"), vmaximum.values),
            "count": (("time"), vcount.values),
            "variance": (("time"), vvariance.values),
            "median": (("time"), vmedian.values),
            "rms": (("time"), vrms.values),            
        },
        {"time": vmean['time']},
    )
    
    ds.to_netcdf(output_filename)