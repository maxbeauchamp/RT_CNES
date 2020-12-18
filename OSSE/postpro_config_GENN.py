#!/usr/bin/env python

__author__ = "Maxime Beauchamp"
__version__ = "0.1"
__date__ = "2020-12-10"
__email__ = "maxime.beauchamp@imt-atantique.fr"

from graphics_OSSE import *

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

NN_lag     = sys.argv[1]
type_obs   = sys.argv[2]
domain     = sys.argv[3] 

workpath = "/users/local/m19beauc/4DVARNN-DinAE_xp/"+domain+"/OSSE/scores_GENN_NNnadlag_"+NN_lag+"_"+type_obs
scratchpath = "/users/local/m19beauc/4DVARNN-DinAE_xp/"+domain+"/OSSE"
if not os.path.exists(workpath):
    mk_dir_recursive(workpath)
#else:
#    shutil.rmtree(workpath)
#    mk_dir_recursive(workpath)    

## parameters
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
#lon = lon[:indLon]
#lat = lat[:indLat]

## store all data in a list
# nadir
FP_GENN_nadir_file_sup1        = scratchpath+'/resIA_nadir_nadlag_'+NN_lag+"_"+type_obs+'/FP_GENN_womissing_wocov/saved_path_000_GENN_womissing.pickle'
FP_GENN_nadir_file_sup2        = scratchpath+'/resIA_nadir_nadlag_'+NN_lag+"_"+type_obs+'/FP_GENN_wmissing_wocov/saved_path_000_GENN_wmissing.pickle'
FP_GENN_nadir_file_unsup       = scratchpath+'/resIA_nadir_nadlag_'+NN_lag+"_"+type_obs+'/FP_GENN_wwmissing_wocov/saved_path_000_GENN_wwmissing.pickle'
FP_GENN_nadir_file_sup1_wOI    = scratchpath+'/resIA_nadir_nadlag_'+NN_lag+"_"+type_obs+'/FP_GENN_womissing_wOI/saved_path_000_GENN_womissing.pickle'
FP_GENN_nadir_file_sup2_wOI    = scratchpath+'/resIA_nadir_nadlag_'+NN_lag+"_"+type_obs+'/FP_GENN_wmissing_wOI/saved_path_000_GENN_wmissing.pickle'
FP_GENN_nadir_file_unsup_wOI   = scratchpath+'/resIA_nadir_nadlag_'+NN_lag+"_"+type_obs+'/FP_GENN_wwmissing_wOI/saved_path_000_GENN_wwmissing.pickle'
# nadir+SWOT
FP_GENN_nadirswot_file_sup1        = scratchpath+'/resIA_nadirswot_nadlag_'+NN_lag+"_"+type_obs+'/FP_GENN_womissing_wocov/saved_path_000_GENN_womissing.pickle'
FP_GENN_nadirswot_file_sup2        = scratchpath+'/resIA_nadirswot_nadlag_'+NN_lag+"_"+type_obs+'/FP_GENN_wmissing_wocov/saved_path_000_GENN_wmissing.pickle'
FP_GENN_nadirswot_file_unsup       = scratchpath+'/resIA_nadirswot_nadlag_'+NN_lag+"_"+type_obs+'/FP_GENN_wwmissing_wocov/saved_path_000_GENN_wwmissing.pickle'
FP_GENN_nadirswot_file_sup1_wOI    = scratchpath+'/resIA_nadirswot_nadlag_'+NN_lag+"_"+type_obs+'/FP_GENN_womissing_wOI/saved_path_000_GENN_womissing.pickle'
FP_GENN_nadirswot_file_sup2_wOI    = scratchpath+'/resIA_nadirswot_nadlag_'+NN_lag+"_"+type_obs+'/FP_GENN_wmissing_wOI/saved_path_000_GENN_wmissing.pickle'
FP_GENN_nadirswot_file_unsup_wOI   = scratchpath+'/resIA_nadirswot_nadlag_'+NN_lag+"_"+type_obs+'/FP_GENN_wwmissing_wOI/saved_path_000_GENN_wwmissing.pickle'

# Reload GENN results
with open(FP_GENN_nadir_file_sup1, 'rb') as handle:
    GT, Obs_nadir, itrp_FP_GENN_nadir_sup1 = pickle.load(handle)[5:8]
with open(FP_GENN_nadir_file_sup2, 'rb') as handle:
    itrp_FP_GENN_nadir_sup2 = pickle.load(handle)[7]
with open(FP_GENN_nadir_file_unsup, 'rb') as handle:
    itrp_FP_GENN_nadir_unsup = pickle.load(handle)[7]
with open(FP_GENN_nadir_file_sup1_wOI, 'rb') as handle:
    itrp_FP_GENN_nadir_sup1_wOI = pickle.load(handle)[7]
with open(FP_GENN_nadir_file_sup2_wOI, 'rb') as handle:
    itrp_FP_GENN_nadir_sup2_wOI = pickle.load(handle)[7]
with open(FP_GENN_nadir_file_sup2_wOI, 'rb') as handle:
    itrp_FP_GENN_nadir_unsup_WOI = pickle.load(handle)[7]
with open(FP_GENN_nadirswot_file_sup1, 'rb') as handle:
    GT, Obs_nadirswot, itrp_FP_GENN_nadirswot_sup1 = pickle.load(handle)[5:8]
with open(FP_GENN_nadirswot_file_sup2, 'rb') as handle:
    itrp_FP_GENN_nadirswot_sup2 = pickle.load(handle)[7]
with open(FP_GENN_nadirswot_file_unsup, 'rb') as handle:
    itrp_FP_GENN_nadirswot_unsup = pickle.load(handle)[7]
with open(FP_GENN_nadirswot_file_sup1_wOI, 'rb') as handle:
    itrp_FP_GENN_nadirswot_sup1_wOI = pickle.load(handle)[7]
with open(FP_GENN_nadirswot_file_sup2_wOI, 'rb') as handle:
    itrp_FP_GENN_nadirswot_sup2_wOI = pickle.load(handle)[7]
with open(FP_GENN_nadirswot_file_sup2_wOI, 'rb') as handle:
    itrp_FP_GENN_nadirswot_unsup_WOI = pickle.load(handle)[7]

## list of dates
lday1 = [ datetime.strftime(datetime.strptime("2012-10-01",'%Y-%m-%d')\
                          + timedelta(days=60+i),"%Y-%m-%d") for i in range(20) ]
lday2 = [ datetime.strftime(datetime.strptime("2012-10-01",'%Y-%m-%d')\
                          + timedelta(days=140+i),"%Y-%m-%d") for i in range(20) ]
lday3 = [ datetime.strftime(datetime.strptime("2012-10-01",'%Y-%m-%d')\
                          + timedelta(days=220+i),"%Y-%m-%d") for i in range(20) ]
lday4 = [ datetime.strftime(datetime.strptime("2012-10-01",'%Y-%m-%d')\
                          + timedelta(days=300+i),"%Y-%m-%d") for i in range(20) ]
lday  = np.concatenate([lday1,lday2,lday3,lday4])
lday2 = [ datetime.strptime(lday[i],'%Y-%m-%d') for i in range(len(lday)) ]

# list_data (nadir)
list_data   = []
list_data.append(GT[:,:indLat,:indLon])
list_data.append(Obs_nadir[:,:indLat,:indLon])
list_data.append(itrp_FP_GENN_nadir_sup1[:,:indLat,:indLon])
list_data.append(itrp_FP_GENN_nadir_sup2[:,:indLat,:indLon])
list_data.append(itrp_FP_GENN_nadir_unsup[:,:indLat,:indLon])
list_data.append(itrp_FP_GENN_nadir_sup1_wOI[:,:indLat,:indLon])
list_data.append(itrp_FP_GENN_nadir_sup2_wOI[:,:indLat,:indLon])
list_data.append(itrp_FP_GENN_nadir_unsup_wOI[:,:indLat,:indLon])

# compare shapes and do appropriate downscaling with minimal resolution
min_res=1e9
for i in range(len(list_data)):
    min_res=min(min_res,list_data[i].shape[1])
for i in range(len(list_data)):
    if list_data[i].shape[1]>min_res:
        dwscale      = int(list_data[i].shape[1]/min_res)
        list_data[i] = einops.reduce(list_data[i], '(t t1) (h h1) (w w1) -> t h w', t1=1, h1=dwscale, w1=dwscale, reduction=np.nanmean)
    print(list_data[i].shape)
dwscale = int(200/min_res)
indLon  = int(indLon/dwscale)
indLat  = int(indLat/dwscale)
lon = np.arange(extent[0],extent[1],1/(20/dwscale))
lat = np.arange(extent[2],extent[3],1/(20/dwscale))

# arguments for plots (nadir)
labels_data = np.array(['GT','Obs','Supervised 1','Supervised 2','Unsupervised','Supervised 1 + OI','Supervised 2 + OI','Unsupervosed + OI'])
colors      = np.array(['k','','red','red','red','blue','blue','blue'])
symbols     = np.array(['k','','o','o','o','o','o','o'])
lstyle      = np.array(['solid','','solid','dashed','dotted','solid','dashed','dotted'])
lwidth      = np.array([2,2,1,1,1,1,1,1])
## Export methods to NetCDF
ncdf_file=workpath+"/NetCDF_nadir_GENNs.nc"
export_NetCDF(list_data,labels_data,lday,lon,lat,ncdf_file)
## nRMSE time series
resfile=workpath+"/TS_nRMSE_nadir_GENNs.png"
plot_nRMSE(list_data,labels_data,colors,symbols,lstyle,lwidth,lday,resfile,gradient=False)

# list_data (nadir+swot)
list_data   = []
list_data.append(GT[:,:indLat,:indLon])
list_data.append(Obs_nadirswot[:,:indLat,:indLon])
list_data.append(itrp_FP_GENN_nadirswot_sup1[:,:indLat,:indLon])
list_data.append(itrp_FP_GENN_nadirswot_sup2[:,:indLat,:indLon])
list_data.append(itrp_FP_GENN_nadirswot_unsup[:,:indLat,:indLon])
list_data.append(itrp_FP_GENN_nadirswot_sup1_wOI[:,:indLat,:indLon])
list_data.append(itrp_FP_GENN_nadirswot_sup2_wOI[:,:indLat,:indLon])
list_data.append(itrp_FP_GENN_nadirswot_unsup_wOI[:,:indLat,:indLon])

# compare shapes and do appropriate downscaling with minimal resolution
min_res=1e9
for i in range(len(list_data)):
    min_res=min(min_res,list_data[i].shape[1])
for i in range(len(list_data)):
    if list_data[i].shape[1]>min_res:
        dwscale      = int(list_data[i].shape[1]/min_res)
        list_data[i] = einops.reduce(list_data[i], '(t t1) (h h1) (w w1) -> t h w', t1=1, h1=dwscale, w1=dwscale, reduction=np.nanmean)
    print(list_data[i].shape)
dwscale = int(200/min_res)
indLon  = int(indLon/dwscale)
indLat  = int(indLat/dwscale)
lon = np.arange(extent[0],extent[1],1/(20/dwscale))
lat = np.arange(extent[2],extent[3],1/(20/dwscale))

## Export methods to NetCDF
ncdf_file=workpath+"/NetCDF_nadirswot_GENNs.nc"
export_NetCDF(list_data,labels_data,lday,lon,lat,ncdf_file)
## nRMSE time series
resfile=workpath+"/TS_nRMSE_nadirswot_GENNs.png"
plot_nRMSE(list_data,labels_data,colors,symbols,lstyle,lwidth,lday,resfile,gradient=False)


