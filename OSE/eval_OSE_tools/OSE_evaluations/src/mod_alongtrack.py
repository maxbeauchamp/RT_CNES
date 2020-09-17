import xarray as xr
import numpy as np
import pandas as pd

def duacs_alongtrack_dataset(files,
                     lon_min=0., lon_max=360., 
                     lat_min=-90, lat_max=90., 
                     time_min='1900-10-01', time_max='2100-01-01'):
    
    def preprocess_duacs_alongtrack(ds_in):
        deltat = ds_in['DeltaT'].values
        dataindexes = ds_in['DataIndexes'].values
    
        time = []
        lon = []
        lat = []
        sla = []
        for cycles in range(len(ds_in.Cycles)):
            begindates = np.ma.masked_invalid(ds_in['BeginDates'].values[:, cycles]).compressed()
            nbpoints = np.ma.masked_invalid(ds_in['NbPoints'][:].values).compressed().astype(np.int64)
            if begindates.size > 0:
                time_alongtrack = np.repeat(begindates, nbpoints, axis=0) + pd.to_timedelta((dataindexes * deltat), unit='s')
                time.append(time_alongtrack)
                lon.append(ds_in['Longitudes'].values)
                lat.append(ds_in['Latitudes'].values)
                sla.append(ds_in['SLA'].values[:, cycles])
        
        time = np.asarray(time).flatten()
        lon = np.asarray(lon).flatten()
        lat = np.asarray(lat).flatten()
        sla = np.asarray(sla).flatten()

        ds_new = xr.Dataset({'SLA': (['time'], sla)},
                            coords={'longitude': (['time'], lon),
                                    'latitude': (['time'], lat),
                                    'time': time})
        
        ds_new = ds_new.sortby(ds_new['time'])
    
        return ds_new

    ds = xr.open_mfdataset(files, drop_variables='Cycles', 
                           preprocess=preprocess_duacs_alongtrack, combine='nested', concat_dim='time', parallel=False)
    
    ds = ds.sel(time=slice(time_min, time_max))
    ds = ds.where((ds["longitude"] >= lon_min) & (ds["longitude"] <= lon_max), drop=True)
    ds = ds.where((ds["latitude"] >= lat_min) & (ds["latitude"] <= lat_max), drop=True)
    
    return ds


def merged_alongtrack(file,
                      lon_min=0., lon_max=360., 
                      lat_min=-90, lat_max=90., 
                      time_min='1900-10-01', time_max='2100-01-01'):
    
    ds = xr.open_dataset(file)
    ds = ds.sel(time=slice(time_min, time_max))
    ds = ds.where((ds["longitude"] >= lon_min) & (ds["longitude"] <= lon_max), drop=True)
    ds = ds.where((ds["latitude"] >= lat_min) & (ds["latitude"] <= lat_max), drop=True)
    
    return ds


def cmems_alongtrack_dataset(dataset,
                      lon_min=0., lon_max=360., 
                      lat_min=-90, lat_max=90., 
                      time_min='1900-10-01', time_max='2100-01-01'):
    
    ds = xr.open_zarr(dataset)
    ds = ds.sel(time=slice(time_min, time_max))
    ds = ds.where((ds["longitude"] >= lon_min) & (ds["longitude"] <= lon_max), drop=True)
    ds = ds.where((ds["latitude"] >= lat_min) & (ds["latitude"] <= lat_max), drop=True)
    ds = ds.rename({'sla_unfiltered': 'SLA'})
    
    return ds

def cmems_alongtrack_dataset_v2(files,
                      lon_min=0., lon_max=360., 
                      lat_min=-90, lat_max=90., 
                      time_min='1900-10-01', time_max='2100-01-01'):
    
    def preprocess_cmems_alongtrack(ds_in):
        ds_new = xr.Dataset({'SLA': (['time'], ds_in.sla_unfiltered.values),
                             'mdt': (['time'], ds_in.mdt.values)},
                             coords={'longitude': (['time'], ds_in.longitude.values),
                                    'latitude': (['time'], ds_in.latitude.values),
                                    'time': ds_in.time.values})
        ds_new = ds_new.sortby(ds_new['time'])
        return ds_new

    ds = xr.open_mfdataset(files,preprocess=preprocess_cmems_alongtrack,\
                            combine='nested', concat_dim='time', parallel=False)    
    ds = ds.sel(time=slice(time_min, time_max))
    ds = ds.where((ds["longitude"] >= lon_min) & (ds["longitude"] <= lon_max), drop=True)
    ds = ds.where((ds["latitude"] >= lat_min) & (ds["latitude"] <= lat_max), drop=True)
    
    return ds
