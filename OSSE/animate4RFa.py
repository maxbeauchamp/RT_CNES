#!/usr/bin/env python

from pb_anda import *
from netCDF4 import Dataset
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

type_obs = sys.argv[1]
domain   = sys.argv[2] 
workpath    = "/gpfsscratch/rech/yrf/uba22to/DINAE_keras/4Ronan/"+domain
scratchpath = '/gpfsscratch/rech/yrf/uba22to/DINAE/'+domain
datapath    = '/gpfswork/rech/yrf/uba22to/DATA/'
if not os.path.exists(workpath):
    mk_dir_recursive(workpath)
else:
    shutil.rmtree(workpath)
    mk_dir_recursive(workpath)

# Reload saved ConvAE and GE-NN results
file_results_nadirswot_Tt=scratchpath+'/resIA_nadirswot_nadlag_5_'+type_obs+'/FP_GENN_wwmissing_wOI/saved_path_019_FP_GENN_wwmissing.pickle'
file_results_nadirswot_Tr=scratchpath+'/resIA_nadirswot_nadlag_5_'+type_obs+'/FP_GENN_wwmissing_wOI/saved_path_019_FP_GENN_wwmissing_train.pickle'
with open(file_results_nadirswot_Tt, 'rb') as handle:
    GT_Tt, Obs_Tt, itrp_FP_GENN_Tt, rec_FP_GENN_Tt, itrp_OI_Tt = pickle.load(handle)
with open(file_results_nadirswot_Tr, 'rb') as handle:
    GT_Tr, Obs_Tr, itrp_FP_GENN_Tr, rec_FP_GENN_Tr, itrp_OI_Tr = pickle.load(handle)

# merge train and test datasets
itrp_OI=np.concatenate((itrp_OI_Tr[0:55],itrp_OI_Tt[0:20],itrp_OI_Tr[55:115],itrp_OI_Tt[20:40],itrp_OI_Tr[115:175],itrp_OI_Tt[40:60],itrp_OI_Tr[175:235],itrp_OI_Tt[60:],itrp_OI_Tr[235:]))
itrp_FP_GENN=np.concatenate((itrp_FP_GENN_Tr[0:55],itrp_FP_GENN_Tt[0:20],itrp_FP_GENN_Tr[55:115],itrp_FP_GENN_Tt[20:40],itrp_FP_GENN_Tr[115:175],itrp_FP_GENN_Tt[40:60],itrp_FP_GENN_Tr[175:235],itrp_FP_GENN_Tt[60:],itrp_FP_GENN_Tr[235:]))

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

# reload GT and obs from original datasets
fileMod=datapath+domain+"/ref/NATL60-CJM165_"+domain+"_ssh_y2013.1y.nc"
fileObs=datapath+domain+"/data/gridded_data_swot_wocorr/dataset_nadir_5d_swot.nc"
nc_data_mod = Dataset(fileMod,'r')
nc_data_obs = Dataset(fileObs,'r')  
GT  = np.copy(nc_data_mod['ssh'][5:360,:,:])
Obs = np.copy(nc_data_obs['ssh_'+type_obs][5:360,:,:])
nc_data_mod.close()
nc_data_obs.close()

print(GT.shape)
# Compute GT gradient
Grad_gt=np.zeros((len(GT),indLat,indLon))
for i in range(len(GT)):
    Grad_gt[i,:,:] = Gradient(GT[i,:indLat,:indLon],2)

def animate(i):
    # Load data
    print(i)
    gt                          = GT[i,:indLat,:indLon]
    Grad_gt                     = Gradient(gt,2)
    OI                          = itrp_OI[i,:indLat,:indLon]
    Grad_OI                     = Gradient(OI,2)
    obs                         = Obs[i,:indLat,:indLon]
    Grad_obs                    = Gradient(obs,2)
    FP_GENN                     = itrp_FP_GENN[i,:indLat,:indLon]
    Grad_FP_GENN                = Gradient(FP_GENN,2)
    for ivar in range(len(var)):   
        ii = int(np.floor(ivar/2))
        jj = int(np.floor(ivar%2))
        ax[ii][jj].clear()
        plot(ax,ii,jj,lon,lat,eval(var[ivar]),title[ivar],\
              extent=extent,cmap=cmap,vmin=vmin,vmax=vmax,colorbar=False)
    return ax[0][0], ax[0][1], ax[1][0], ax[1][1]

# for nadirswot
var=['gt','obs','OI','FP_GENN']
title=['Ground Truth','Obs','OI','FP-GENN']
vmax=np.round(np.nanmax([np.abs(GT),np.abs(GT)]),2)
vmin=-1.0*vmax
cmap="coolwarm"
fig, ax = plt.subplots(2,2,figsize=(15,15),\
          subplot_kw=dict(projection=ccrs.PlateCarree(central_longitude=0.0)))
plt.subplots_adjust(hspace=0.5)
ani = animation.FuncAnimation(fig, animate, frames=np.arange(1,355), blit=False, interval=1000, repeat=False)
writer = animation.FFMpegWriter(fps=3, bitrate=5000)
ani.save(workpath+"/animation_maps.mp4", writer = writer)
plt.close()

var=['Grad_gt','Grad_obs','Grad_OI','Grad_FP_GENN']
title=[r"$\nabla_{GT}$",r"$\nabla_{Obs}$",r"$\nabla_{OI}$",r"$\nabla_{FP-GENN}$"]
vmax=np.round(np.nanmax([np.abs(Grad_gt),np.abs(Grad_gt)]),2)
vmin=0.
cmap="viridis"
fig, ax = plt.subplots(2,2,figsize=(15,15),\
          subplot_kw=dict(projection=ccrs.PlateCarree(central_longitude=0.0)))
plt.subplots_adjust(hspace=0.5)
ani = animation.FuncAnimation(fig, animate, frames=np.arange(1,355), blit=False, interval=1000, repeat=False)
writer = animation.FFMpegWriter(fps=3, bitrate=5000)
ani.save(workpath+"/animation_grads.mp4", writer = writer)
plt.close()


