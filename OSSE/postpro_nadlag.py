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

type_obs   = sys.argv[1]
domain     = sys.argv[2] 

workpath = "/users/local/m19beauc/4DVARNN-DinAE_xp/"+domain+"/OSSE/scores_allmethods_nadlag_"+type_obs
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
AnDA_nadir_lag_0_file             = scratchpath+'/resAnDA_nadir_nadlag_0_'+type_obs+'/saved_path.pickle'
FP_GENN_nadir_lag_0_file          = scratchpath+'/resIA_nadir_nadlag_0'+type_obs+'/FP_GENN_wmissing_wOI/saved_path_000_GENN_wmissing.pickle'
AnDA_nadir_lag_5_file             = scratchpath+'/resAnDA_nadir_nadlag_5_'+type_obs+'/saved_path.pickle'
FP_GENN_nadir_lag_5_file          = scratchpath+'/resIA_nadir_nadlag_5'+type_obs+'/FP_GENN_wmissing_wOI/saved_path_000_GENN_wmissing.pickle'
AnDA_nadirswot_lag_0_file         = scratchpath+'/resAnDA_nadirswot_nadlag_0'+type_obs+'/saved_path.pickle'
FP_GENN_nadirswot_lag_0_file      = scratchpath+'/resIA_nadirswot_nadlag_0'+type_obs+'/FP_GENN_wmissing_wOI/saved_path_000_GENN_wmissing.pickle'
AnDA_nadirswot_lag_5_file         = scratchpath+'/resAnDA_nadirswot_nadlag_5'+type_obs+'/saved_path.pickle'
FP_GENN_nadirswot_lag_5_file      = scratchpath+'/resIA_nadirswot_nadlag_5'+type_obs+'/FP_GENN_wmissing_wOI/saved_path_000_GENN_wmissing.pickle'

# Reload saved AnDA result
with open(AnDA_nadir_lag_0_file, 'rb') as handle:
    AnDA_ssh_1, itrp_dineof = pickle.load(handle)
    AnDA_ssh_1_nadir_0 = AnDA_ssh_1  
    itrp_dineof_nadir_0 = itrp_dineof
with open(AnDA_nadirswot_lag_0_file, 'rb') as handle:
    AnDA_ssh_1, itrp_dineof = pickle.load(handle)
    AnDA_ssh_1_nadirswot_0 = AnDA_ssh_1      
    itrp_dineof_nadirswot_0 = itrp_dineof
with open(AnDA_nadir_lag_5_file, 'rb') as handle:
    AnDA_ssh_1, itrp_dineof = pickle.load(handle)
    AnDA_ssh_1_nadir_5 = AnDA_ssh_1
    itrp_dineof_nadir_5 = itrp_dineof
with open(AnDA_nadirswot_lag_5_file, 'rb') as handle:
    AnDA_ssh_1, itrp_dineof = pickle.load(handle)
    AnDA_ssh_1_nadirswot_5 = AnDA_ssh_1
    itrp_dineof_nadirswot_5 = itrp_dineof
# Reload saved ConvAE and GE-NN results
with open(FP_GENN_nadir_lag_0_file, 'rb') as handle:
    itrp_FP_GENN_nadir_0, rec_FP_GENN_nadir_0 = pickle.load(handle)[7:9]
with open(FP_GENN_nadirswot_lag_0_file, 'rb') as handle:
    itrp_FP_GENN_nadirswot_0, rec_FP_GENN_nadirswot_0 = pickle.load(handle)[7:9]
with open(FP_GENN_nadir_lag_5_file, 'rb') as handle:
    itrp_FP_GENN_nadir_5, rec_FP_GENN_nadir_5 = pickle.load(handle)[7:9]
with open(FP_GENN_nadirswot_lag_5_file, 'rb') as handle:
    itrp_FP_GENN_nadirswot_5, rec_FP_GENN_nadirswot_5 = pickle.load(handle)[7:9]


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

GT = AnDA_ssh_1_nadir.GT[:,:indLat,:indLon]
# list_data (AnDA nadir)
list_data   = []
list_data.append(GT)
list_data.append(AnDA_ssh_1_nadir_0.itrp_postAnDA[:,:indLat,:indLon])
list_data.append(AnDA_ssh_1_nadir_5.itrp_postAnDA[:,:indLat,:indLon])
# arguments for plots (nadir)
labels_data = np.array(['GT','Obs','Post-AnDA (lag=0)','Post-AnDA (lag=5)'])
colors      = np.array(['k','','red','blue'])
symbols     = np.array(['k','','o','o'])
lstyle      = np.array(['solid','','solid','solid'])
lwidth      = np.array([2,2,1,1])
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
## nRMSE time series
resfile=workpath+"/TS_AnDA_nadir_nadlag.png"
plot_nRMSE(list_data,labels_data,colors,symbols,lstyle,lwidth,lday,resfile,gradient=False)

# list_data (AnDA nadirswot)
list_data   = []
list_data.append(GT)
list_data.append(AnDA_ssh_1_nadirswot_0.itrp_postAnDA[:,:indLat,:indLon])
list_data.append(AnDA_ssh_1_nadirswot_5.itrp_postAnDA[:,:indLat,:indLon])
# arguments for plots (nadirswot)
labels_data = np.array(['GT','Obs','Post-AnDA (lag=0)','Post-AnDA (lag=5)'])
colors      = np.array(['k','','red','blue'])
symbols     = np.array(['k','','o','o'])
lstyle      = np.array(['solid','','solid','solid'])
lwidth      = np.array([2,2,1,1])
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
## nRMSE time series
resfile=workpath+"/TS_AnDA_nadirswot_nadlag.png"
plot_nRMSE(list_data,labels_data,colors,symbols,lstyle,lwidth,lday,resfile,gradient=False)

# list_data (GENN nadir)
list_data   = []
list_data.append(GT)
list_data.append(itrp_FP_GENN_nadir_0[:,:indLat,:indLon])
list_data.append(itrp_FP_GENN_nadir_5[:,:indLat,:indLon])
# arguments for plots (nadir)
labels_data = np.array(['GT','Obs','FP-GENN (lag=0)','FP-GENN (lag=5)'])
colors      = np.array(['k','','red','blue'])
symbols     = np.array(['k','','o','o'])
lstyle      = np.array(['solid','','solid','solid'])
lwidth      = np.array([2,2,1,1])
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
## nRMSE time series
resfile=workpath+"/TS_GENN_nadir_nadlag.png"
plot_nRMSE(list_data,labels_data,colors,symbols,lstyle,lwidth,lday,resfile,gradient=False)

# list_data (GENN nadirswot)
list_data   = []
list_data.append(GT)
list_data.append(itrp_FP_GENN_nadirswot_0[:,:indLat,:indLon])
list_data.append(itrp_FP_GENN_nadirswot_5[:,:indLat,:indLon])
# arguments for plots (nadirswot)
labels_data = np.array(['GT','Obs','FP-GENN (lag=0)','FP-GENN (lag=5)'])
colors      = np.array(['k','','red','blue'])
symbols     = np.array(['k','','o','o'])
lstyle      = np.array(['solid','','solid','solid'])
lwidth      = np.array([2,2,1,1])
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
## nRMSE time series
resfile=workpath+"/TS_GENN_nadirswot_nadlag.png"
plot_nRMSE(list_data,labels_data,colors,symbols,lstyle,lwidth,lday,resfile,gradient=False)

