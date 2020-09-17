import xarray as xr
import numpy as np
import pyinterp
import netCDF4


def duacs_grid_dataset(list_of_file, variable_name='Grid_0001', 
                       lon_min=0., lon_max=360., 
                       lat_min=-90, lat_max=90., 
                       time_min='1900-10-01', time_max='2100-01-01', is_circle=True):
    
    def preprocess_duacs_maps(ds):
        
        vtime = ds[variable_name].attrs['Date_CNES_JD']        
        ds.coords['time'] = np.datetime64(netCDF4.num2date(vtime, units='days since 1950-01-01'))
        
        if ds[variable_name].units == 'cm':
            ds[variable_name] = ds[variable_name] / 100.
            ds[variable_name].attrs = {'units': 'm'}

        return ds
    
    ds = xr.open_mfdataset(list_of_file, concat_dim ='time', combine='nested', parallel=True, preprocess=preprocess_duacs_maps)
    ds = ds.sel(time=slice(time_min, time_max))
    ds = ds.where((ds["NbLongitudes"] >= lon_min) & (ds["NbLongitudes"] <= lon_max), drop=True)
    ds = ds.where((ds["NbLatitudes"] >= lat_min) & (ds["NbLatitudes"] <= lat_max), drop=True)
    
    x_axis = pyinterp.Axis(ds["NbLongitudes"][:], is_circle=is_circle)
    y_axis = pyinterp.Axis(ds["NbLatitudes"][:])
    z_axis = pyinterp.TemporalAxis(ds["time"][:])
    
    var = ds[variable_name][:].transpose('NbLongitudes', 'NbLatitudes', 'time')
    # The undefined values must be set to nan.
    try:
        var[var.mask] = float("nan")
    except AttributeError:
        pass
    
    grid = pyinterp.Grid3D(x_axis, y_axis, z_axis, var.data)
    
    del ds
    
    return x_axis, y_axis, z_axis, grid


def miost_grid_dataset(list_of_file, var2add, var2sub, 
                       lon_min=0., lon_max=360., 
                       lat_min=-90, lat_max=90., 
                       time_min='1900-10-01', time_max='2100-01-01', is_circle=True):
    
    
    ds = xr.open_mfdataset(list_of_file, concat_dim ='time', combine='nested', parallel=True)
    ds = ds.sel(time=slice(time_min, time_max))
    ds = ds.where((ds["longitude"] >= lon_min) & (ds["longitude"] <= lon_max), drop=True)
    ds = ds.where((ds["latitude"] >= lat_min) & (ds["latitude"] <= lat_max), drop=True)
    
    # print(ds)
    
    x_axis = pyinterp.Axis(ds["longitude"][:], is_circle=is_circle)
    y_axis = pyinterp.Axis(ds["latitude"][:])
    z_axis = pyinterp.TemporalAxis(ds["time"][:])
    
    for variable_name in var2add:
        try:
            var += ds[variable_name][:]   #ds['Hb'][:] + ds['Hls'][:] + ds['Hss'][:] # + ds['Heq'][:]
        except UnboundLocalError:
            var = ds[variable_name][:]
    
    for variable_name in var2sub:
        try:
            var -= ds[variable_name][:] 
        except UnboundLocalError:
            var = ds[variable_name][:]
        
    var = var.transpose('longitude', 'latitude', 'time')
    # The undefined values must be set to nan.
    try:
        var[var.mask] = float("nan")
    except AttributeError:
        pass
    
    grid = pyinterp.Grid3D(x_axis, y_axis, z_axis, var.data)
    
    del ds
    
    return x_axis, y_axis, z_axis, grid    


def dymost_grid_dataset(list_of_file, var2add, var2sub, 
                       lon_min=0., lon_max=360., 
                       lat_min=-90, lat_max=90., 
                       time_min='1900-10-01', time_max='2100-01-01', is_circle=True):
    
    
    ds = xr.open_mfdataset(list_of_file, concat_dim ='time', combine='nested', parallel=True)
    ds = ds.sel(time=slice(time_min, time_max))
    ds = ds.where((ds["lon"] >= lon_min) & (ds["lon"] <= lon_max), drop=True)
    ds = ds.where((ds["lat"] >= lat_min) & (ds["lat"] <= lat_max), drop=True)
    
    # print(ds)
    
    x_axis = pyinterp.Axis(ds["lon"][0, :], is_circle=is_circle)
    y_axis = pyinterp.Axis(ds["lat"][:, 0])
    z_axis = pyinterp.TemporalAxis(ds["time"][:])
    
    for variable_name in var2add:
        try:
            var += ds[variable_name][:]   #ds['Ha'][:]
        except UnboundLocalError:
            var = ds[variable_name][:]
    
    for variable_name in var2sub:
        try:
            var -= ds[variable_name][:] 
        except UnboundLocalError:
            var = ds[variable_name][:]
        
    var = var.transpose('x', 'y', 'time')
    # The undefined values must be set to nan.
    try:
        var[var.mask] = float("nan")
    except AttributeError:
        pass
    
    grid = pyinterp.Grid3D(x_axis, y_axis, z_axis, var.data)
    
    del ds
    
    return x_axis, y_axis, z_axis, grid 


def fpgenn_grid_dataset(list_of_file, var2add, var2sub, 
                       lon_min=0., lon_max=360., 
                       lat_min=-90, lat_max=90., 
                       time_min='1900-10-01', time_max='2100-01-01', is_circle=True):
    
    
    ds = xr.open_mfdataset(list_of_file, concat_dim ='time', combine='nested', parallel=True)
    ds = ds.sel(time=slice(time_min, time_max))
    ds = ds.where((ds["lon"]%360. >= lon_min) & (ds["lon"]%360. <= lon_max), drop=True)
    ds = ds.where((ds["lat"] >= lat_min) & (ds["lat"] <= lat_max), drop=True)
    
    # print(ds)
    
    x_axis = pyinterp.Axis(ds["lon"][:]%360., is_circle=is_circle)
    y_axis = pyinterp.Axis(ds["lat"][:])
    z_axis = pyinterp.TemporalAxis(ds["time"][:])
    
    for variable_name in var2add:
        try:
            var += ds[variable_name][:]   #ds['Ha'][:]
        except UnboundLocalError:
            var = ds[variable_name][:]
    
    for variable_name in var2sub:
        try:
            var -= ds[variable_name][:] 
        except UnboundLocalError:
            var = ds[variable_name][:]

    # MB clean boundary for file OSE_GULFSTREAM_FPGENN.nc
    #var.values[:, 0:3, :] = np.nan
    #var.values[:, :, 0:3] = np.nan
    
    var = var.transpose('lon', 'lat', 'time')
    # The undefined values must be set to nan.
    try:
        var[var.mask] = float("nan")
    except AttributeError:
        pass
    
    grid = pyinterp.Grid3D(x_axis, y_axis, z_axis, var.data)
    
    del ds
    
    return x_axis, y_axis, z_axis, grid


def bfn_grid_dataset(list_of_file, var2add, var2sub, 
                       lon_min=0., lon_max=360., 
                       lat_min=-90, lat_max=90., 
                       time_min='1900-10-01', time_max='2100-01-01', is_circle=True):
    
    
    ds = xr.open_mfdataset(list_of_file, concat_dim ='time', combine='nested', parallel=True)
    ds = ds.sel(time=slice(time_min, time_max))
    ds = ds.where((ds["lon"]%360. >= lon_min) & (ds["lon"]%360. <= lon_max), drop=True)
    ds = ds.where((ds["lat"] >= lat_min) & (ds["lat"] <= lat_max), drop=True)
    
    x_axis = pyinterp.Axis(ds["lon"][:]%360., is_circle=is_circle)
    y_axis = pyinterp.Axis(ds["lat"][:])
    z_axis = pyinterp.TemporalAxis(ds["time"][:])
    
    for variable_name in var2add:
        try:
            var += ds[variable_name][:]
        except UnboundLocalError:
            var = ds[variable_name][:]
    
    for variable_name in var2sub:
        try:
            var -= ds[variable_name][:] 
        except UnboundLocalError:
            var = ds[variable_name][:]

    # MB clean boundary for file OSE_GULFSTREAM_FPGENN.nc
    #var.values[:, 0:3, :] = np.nan
    #var.values[:, :, 0:3] = np.nan
    
    # ds['time'] = (ds['time'] - np.datetime64('1950-01-01T00:00:00Z')) / np.timedelta64(1, 'D')
    
    var = var.transpose('lon', 'lat', 'time')

    # The undefined values must be set to nan.
    try:
        var[var.mask] = float("nan")
    except AttributeError:
        pass
    
    grid = pyinterp.Grid3D(x_axis, y_axis, z_axis, var.data)
    
    del ds
    
    return x_axis, y_axis, z_axis, grid