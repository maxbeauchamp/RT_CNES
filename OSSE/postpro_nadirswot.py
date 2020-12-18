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

AnDA_lag   = sys.argv[1]
NN_lag     = sys.argv[2]
type_obs   = sys.argv[3]
domain     = sys.argv[4] 

workpath = "/users/local/m19beauc/4DVARNN-DinAE_xp/"+domain+"/OSSE/scores_allmethods_AnDAnadlag_"+AnDA_lag+"_NNnadlag_"+NN_lag+"_"+type_obs
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
AnDA_nadirswot_file             = scratchpath+'/resAnDA_nadirswot_nadlag_'+AnDA_lag+"_"+type_obs+'/saved_path.pickle'
FP_ConvAE_nadirswot_file        = scratchpath+'/resIA_nadirswot_nadlag_'+NN_lag+"_"+type_obs+'/FP_ConvAE_wmissing_wOI/saved_path_000_ConvAE_wmissing.pickle'
FP_GENN_nadirswot_file_sup      = scratchpath+'/resIA_nadirswot_nadlag_'+NN_lag+"_"+type_obs+'/FP_GENN_wmissing_wOI/saved_path_000_GENN_wmissing.pickle'
FP_GENN_nadirswot_file_unsup    = scratchpath+'/resIA_nadirswot_nadlag_'+NN_lag+"_"+type_obs+'/FP_GENN_wwmissing_wOI/saved_path_019_GENN_wwmissing.pickle'
Grad_GENN_nadirswot_file_sup    = scratchpath+'/resIA_nadirswot_nadlag_'+NN_lag+"_"+type_obs+'/GB1_GENN_wmissing_wOI/saved_path_001_GENN_wmissing.pickle'
Grad_GENN_nadirswot_file_unsup  = scratchpath+'/resIA_nadirswot_nadlag_'+NN_lag+"_"+type_obs+'/GB1_GENN_wwmissing_wOI/saved_path_000_GENN_wwmissing.pickle'

# Reload saved AnDA result
with open(AnDA_nadirswot_file, 'rb') as handle:
    AnDA_ssh_1, itrp_dineof = pickle.load(handle)
    AnDA_ssh_1_nadirswot = AnDA_ssh_1  
    itrp_dineof_nadirswot = itrp_dineof
# Reload saved ConvAE and GE-NN results
with open(FP_ConvAE_nadirswot_file, 'rb') as handle:
    itrp_FP_ConvAE_nadirswot, rec_FP_ConvAE_nadirswot = pickle.load(handle)[7:9]
with open(FP_GENN_nadirswot_file_sup, 'rb') as handle:
    itrp_FP_GENN_nadirswot, rec_FP_GENN_nadirswot = pickle.load(handle)[7:9]
with open(Grad_GENN_nadirswot_file_sup, 'rb') as handle:
    itrp_GB_GENN_nadirswot, rec_GB_GENN_nadirswot = pickle.load(handle)[7:9]

# list_data (nadir+swot)
list_data   = []
list_data.append(AnDA_ssh_1_nadirswot.GT[:,:indLat,:indLon])
list_data.append(AnDA_ssh_1_nadirswot.Obs[:,:indLat,:indLon])
list_data.append(AnDA_ssh_1_nadirswot.itrp_OI[:,:indLat,:indLon])
list_data.append(AnDA_ssh_1_nadirswot.itrp_postAnDA[:,:indLat,:indLon])
list_data.append(itrp_dineof_nadirswot[:,:indLat,:indLon])
#list_data.append(itrp_FP_ConvAE_nadirswot[:,:indLat,:indLon])
list_data.append(itrp_FP_GENN_nadirswot[:,:indLat,:indLon])
list_data.append(itrp_GB_GENN_nadirswot[:,:indLat,:indLon])
# arguments for plots (nadir+swot)
labels_data = np.array(['GT','Obs (nadir+swot)','OI (nadir+swot)','Post-AnDA (nadir+swot)','VE-DINEOF (nadir+swot)','FP-ConvAE (nadir+swot)','FP-GENN (nadir+swot)','Grad-GENN (nadir+swot)'])
list_suffix = np.array(['GT','Obs_nadirswot','OI_nadirswot','Post_AnDA_nadirswot','VE_DINEOF_nadirswot','FP_ConvAE_nadirswot','FP_GENN_nadirswot','Grad_GENN_nadirswot'])
colors      = np.array(['k','','red','seagreen','steelblue','violet','darkorange','darkblue'])
symbols     = np.array(['k','','o','o','o','o','o','o'])
lstyle      = np.array(['solid','','solid','solid','solid','solid','solid','solid'])
lwidth      = np.array([2,2,2,1,1,1,1,1])

labels_data = np.array(['GT','Obs (nadir+swot)','OI (nadir+swot)','Post-AnDA (nadir+swot)','VE-DINEOF (nadir+swot)','FP-GENN (nadir+swot)','Grad-GENN (nadir+swot)'])
list_suffix = np.array(['GT','Obs_nadirswot','OI_nadirswot','Post_AnDA_nadirswot','VE_DINEOF_nadirswot','FP_GENN_nadirswot','Grad_GENN_nadirswot'])
colors      = np.array(['k','','red','seagreen','steelblue','darkorange','darkblue'])
symbols     = np.array(['k','','o','o','o','o','o'])
lstyle      = np.array(['solid','','solid','solid','solid','solid','solid'])
lwidth      = np.array([2,2,2,1,1,1,1])

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

## Export methods to NetCDF
ncdf_file=workpath+"/NetCDF_nadirswot.nc"
export_NetCDF(list_data,labels_data,lday,lon,lat,ncdf_file)
## test PSD
resfile=workpath+"/BOOST_PSD_nadirswot"
plot_psd(ncdf_file,labels_data,lday,resfile)
## Compute R/I/AE scores
resfile=workpath+"/RIAE_scores_nadirswot.png"
RIAE_scores(list_data,labels_data,resfile,gradient=False)
## nRMSE time series
resfile=workpath+"/TS_nRMSE_nadirswot.png"
plot_nRMSE(list_data,labels_data,colors,symbols,lstyle,lwidth,lday,resfile,gradient=False)
resfile=workpath+"/TS_nRMSE_Grad_nadirswot.png"
plot_nRMSE(list_data,labels_data,colors,symbols,lstyle,lwidth,lday,resfile,gradient=True)
## average SNR
resfile=workpath+"/SNR_nadirswot.png"
resssh=4*dwscale
plot_SNR(list_data,labels_data,colors,symbols,lstyle,lwidth,lday,resssh,resfile)
## average Taylor diagrams
resfile=workpath+"/Taylor_diagram_nadirswot.png"
Taylor_diagram(list_data,labels_data,colors,symbols,resfile)
## plot individual maps (SSH & Gradients)
plot_maps(list_data,list_suffix,labels_data,lday,"2013-08-04",extent,lon,lat,workpath)

