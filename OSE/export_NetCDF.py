#!/usr/bin/env python

from netCDF4 import Dataset
import numpy as np
import sys
import pickle
import xarray as xr
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.animation as animation
import mpl_toolkits.mplot3d.axes3d as p3
from mpl_toolkits.axes_grid1.inset_locator import InsetPosition
np.random.seed(1)

# function to create recursive paths
def mk_dir_recursive(dir_path):
    if os.path.isdir(dir_path):
        return
    h, t = os.path.split(dir_path)  # head/tail
    if not os.path.isdir(h):
        mk_dir_recursive(h)

    new_path = join_paths(h, t)
    if not os.path.isdir(new_path):
        os.mkdir(new_path)

def str2bool(v):
  return v.lower() in ("yes", "true", "t", "1")

# setup learning configuration
lag        = sys.argv[1]
domain     = sys.argv[2]
wregul     = str2bool(sys.argv[3])
if wregul==True:
    scratchpath = '/users/local/m19beauc/4DVARNN-DinAE_xp/'+domain+'/OSE/resIA_nadir_nadlag_'+lag+'_obs/GB1_GENN_wwmissing_wOI_wtrain_wregul'
else:
    scratchpath = '/users/local/m19beauc/4DVARNN-DinAE_xp/'+domain+'/OSE/resIA_nadir_nadlag_'+lag+'_obs/GB1_GENN_wwmissing_wOI_wtrain'

datapath    = '/users/local/DATA/OSE/'+domain+'/training'

@staticmethod
def preprocess2(ds):
    ds.time.attrs['units'] = 'seconds since 2017-01-01'
    ds.time.attrs['calendar'] = 'standard'
    ds = xr.decode_cf(ds)
    return ds

# Reload saved GE-NN results
file_results=scratchpath+'/saved_path_005_GENN_wwmissing.pickle'
with open(file_results, 'rb') as handle: 
    Obs, FP_GENN, ref_FP_GENN, OI = pickle.load(handle)[1:5]

file_data=datapath+"/dataset_nadir_0d.nc"
ds = xr.open_dataset(file_data)

# 4-plots video individual maps
if domain=="OSMOSIS":
    extent     = [-19.5,-11.5,45.,55.]
    indLat     = 200
    indLon     = 160
elif domain=='GULFSTREAM':
    extent     = [-65.,-55.,33.,43.]
    indLat     = 200
    indLon     = 200
else:
    extent=[-65.,-55.,30.,40.]
    indLat     = 200
    indLon     = 200
dwscale = int(indLon/FP_GENN.shape[2])
indLon=int(indLon/dwscale)
indLat=int(indLat/dwscale)
lon = np.arange(extent[0],extent[1],1/(20/dwscale))
lat = np.arange(extent[2],extent[3],1/(20/dwscale))
lon = lon[:indLon]
lat = lat[:indLat]
mesh_lat, mesh_lon = np.meshgrid(lat, lon)
time_u = ds.time.values

#for i in range(len(FP_GENN)):
#    FP_GENN[i] = median_filter(FP_GENN[i], size=5)

xrdata = xr.Dataset(\
                data_vars={'longitude': (('lat','lon'),mesh_lon),\
                           'latitude' : (('lat','lon'),mesh_lat),\
                           'Time'     : (('time'),time_u),\
                           'obs'      : (('time','lat','lon'),Obs),\
                           'OI'       : (('time','lat','lon'),OI),\
                           'GB-GENN'  : (('time','lat','lon'),FP_GENN)},\
                coords={'lon': lon,'lat': lat,'time': range(0,len(time_u))})
xrdata.time.attrs['units']='days since 2017-01-01 00:00:00'
xrdata.to_netcdf(path=scratchpath+"/OSE_"+domain+"_GB1_GENN.nc", mode='w')
