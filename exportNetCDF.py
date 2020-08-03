#!/usr/bin/env python

from pb_anda import *
from netCDF4 import Dataset
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

domain   = sys.argv[1] 
workpath    = "/gpfsscratch/rech/yrf/uba22to/DINAE_keras/4MBallarotta/"+domain
scratchpath = '/gpfsscratch/rech/yrf/uba22to/DINAE_keras/OSE/'+domain+\
              '/resIA_nadir_nadlag_0_obs/FP_GENN_wwmissing_wOI_wtrain'
datapath    = '/gpfswork/rech/yrf/uba22to/DATA/OSE/'+domain+'/training/data'

@staticmethod
def preprocess2(ds):
    ds.time.attrs['units'] = 'seconds since 2017-01-01'
    ds.time.attrs['calendar'] = 'standard'
    ds = xr.decode_cf(ds)
    return ds


if not os.path.exists(workpath):
    mk_dir_recursive(workpath)
else:
    shutil.rmtree(workpath)
    mk_dir_recursive(workpath)

# Reload saved GE-NN results
file_results=scratchpath+'/saved_path_FP_GENN_wwmissing.pickle'
with open(file_results, 'rb') as handle: 
    OI, Obs, FP_GENN = pickle.load(handle)[1:4]

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
lon = np.arange(extent[0],extent[1],1/20)
lat = np.arange(extent[2],extent[3],1/20)
lon = lon[:indLon]
lat = lat[:indLat]
mesh_lat, mesh_lon = np.meshgrid(lat, lon)
time_u = ds.time.values
xrdata = xr.Dataset(\
                data_vars={'longitude': (('lat','lon'),mesh_lon),\
                           'latitude' : (('lat','lon'),mesh_lat),\
                           'Time'     : (('time'),time_u),\
                           'obs'      : (('time','lat','lon'),Obs),\
                           'OI'       : (('time','lat','lon'),OI),\
                           'FP-GENN'  : (('time','lat','lon'),FP_GENN)},\
                coords={'lon': lon,'lat': lat,'time': range(0,len(time_u))})
xrdata.time.attrs['units']='days since 2017-01-01 00:00:00'
xrdata.to_netcdf(path=workpath+"/OSE_"+domain+"_FPGENN.nc", mode='w')
