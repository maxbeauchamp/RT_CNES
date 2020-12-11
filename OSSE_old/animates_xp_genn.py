#!/usr/bin/env python

""" animate_results_final.py: script to animate maps MS-PB-AnDA and NN-based methods with three different datasets (nadir / swot / nadirswot) """

__author__ = "Maxime Beauchamp"
__version__ = "0.1"
__date__ = "2020-02-01"
__email__ = "maxime.beauchamp76@gmail.com"

from pb_anda import *
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
workpath = "/users/local/m19beauc/4DVARNN-DinAE_xp/OSSE/"+domain+"/animates_GENN_wcovar_"+type_obs
scratchpath = '/users/local/m19beauc/4DVARNN-DinAE_xp/OSSE/'+domain
if not os.path.exists(workpath):
    mk_dir_recursive(workpath)
else:
    shutil.rmtree(workpath)
    mk_dir_recursive(workpath)

# Reload AnDA results (for GT and OI)
file_results_AnDA=scratchpath+'/resAnDA_nadirswot_nadlag_5_'+type_obs+'/saved_path.pickle'
with open(file_results_AnDA, 'rb') as handle:
    AnDA_ssh_1, itrp_dineof = pickle.load(handle)
    AnDA_ssh_1_nadirswot = AnDA_ssh_1
    itrp_dineof_nadirswot = itrp_dineof

# Reload GE-NN results
file_results_1=scratchpath+'/resIA_tests_GENN_'+type_obs+'/FP_GENN_womissing_wocov/saved_path_019_FP_GENN_womissing.pickle'
with open(file_results_1, 'rb') as handle:
    itrp_FP_GENN_1 = pickle.load(handle)[2]
file_results_2=scratchpath+'/resIA_tests_GENN_'+type_obs+'/FP_GENN_wmissing_wocov/saved_path_019_FP_GENN_wmissing.pickle'
with open(file_results_2, 'rb') as handle:
    itrp_FP_GENN_2 = pickle.load(handle)[2]
file_results_3=scratchpath+'/resIA_tests_GENN_'+type_obs+'/FP_GENN_womissing_wOI/saved_path_019_FP_GENN_womissing.pickle'
with open(file_results_3, 'rb') as handle:
    itrp_FP_GENN_3 = pickle.load(handle)[2]
file_results_4=scratchpath+'/resIA_tests_GENN_'+type_obs+'/FP_GENN_wmissing_wOI/saved_path_019_FP_GENN_wmissing.pickle'
with open(file_results_4, 'rb') as handle:
    itrp_FP_GENN_4 = pickle.load(handle)[2]


# 6-plots video individual maps
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

# Compute GT gradient
Grad_gt=np.zeros((len(AnDA_ssh_1.GT),indLon,indLat))
for i in range(len(AnDA_ssh_1.GT)):
    Grad_gt[i,:,:] = Gradient(AnDA_ssh_1.GT[i,:indLat,:indLon],2)

def init():
    global fig, ax
    # Load data
    gt                          = AnDA_ssh_1.GT[0,:indLat,:indLon]
    Grad_gt                     = Gradient(gt,2)
    obs                         = AnDA_ssh_1_nadirswot.Obs[0,:indLat,:indLon]
    Grad_obs                    = Gradient(obs,2)
    # nadirswot
    FP_GENN_1                   = itrp_FP_GENN_1[0,:indLat,:indLon]
    Grad_FP_GENN_1              = Gradient(FP_GENN_1,2)
    FP_GENN_2                   = itrp_FP_GENN_2[0,:indLat,:indLon]
    Grad_FP_GENN_2              = Gradient(FP_GENN_2,2)
    FP_GENN_3                   = itrp_FP_GENN_3[0,:indLat,:indLon]
    Grad_FP_GENN_3              = Gradient(FP_GENN_3,2)
    FP_GENN_4                   = itrp_FP_GENN_4[0,:indLat,:indLon]
    Grad_FP_GENN_4              = Gradient(FP_GENN_4,2)
    for ivar in range(len(var)):
        ii = int(np.floor(ivar/3))
        jj = int(np.floor(ivar%3))
        plot(ax,ii,jj,lon,lat,eval(var[ivar]),title[ivar],\
              extent=extent,cmap=cmap,vmin=vmin,vmax=vmax,colorbar=True,\
              orientation="vertical")

def animate(i):
    global fig, ax
    # Load data
    gt                          = AnDA_ssh_1.GT[i,:indLat,:indLon]
    Grad_gt                     = Gradient(gt,2)
    obs                         = AnDA_ssh_1_nadirswot.Obs[i,:indLat,:indLon]
    Grad_obs                    = Gradient(obs,2)
    # nadirswot
    FP_GENN_1                   = itrp_FP_GENN_1[i,:indLat,:indLon]
    Grad_FP_GENN_1              = Gradient(FP_GENN_1,2)
    FP_GENN_2                   = itrp_FP_GENN_2[i,:indLat,:indLon]
    Grad_FP_GENN_2              = Gradient(FP_GENN_2,2)
    FP_GENN_3                   = itrp_FP_GENN_3[i,:indLat,:indLon]
    Grad_FP_GENN_3              = Gradient(FP_GENN_3,2)
    FP_GENN_4                   = itrp_FP_GENN_4[i,:indLat,:indLon]
    Grad_FP_GENN_4              = Gradient(FP_GENN_4,2)
    for ivar in range(len(var)):   
        ii = int(np.floor(ivar/3))
        jj = int(np.floor(ivar%3))
        ax[ii][jj].cla()
        plot(ax,ii,jj,lon,lat,eval(var[ivar]),title[ivar],\
              extent=extent,cmap=cmap,vmin=vmin,vmax=vmax,colorbar=False)
    return fig, ax

# for nadir
var=['gt','obs',\
     'FP_GENN_1','FP_GENN_2',\
     'FP_GENN_3','FP_GENN_4']
title=['Ground Truth','Obs (nadir)',\
       'GENN (all data)','GENN (masked data)',\
       'GENN (all data + OI)','GENN (masked data + OI)']
vmax=np.round(np.nanmax([np.abs(AnDA_ssh_1.GT),np.abs(AnDA_ssh_1.GT)]),2)
vmin=-1.0*vmax
cmap="coolwarm"
fig, ax = plt.subplots(2,3,figsize=(15,10),\
          subplot_kw=dict(projection=ccrs.PlateCarree(central_longitude=0.0)))
plt.subplots_adjust(hspace=0.5)
ani = animation.FuncAnimation(fig, animate, init_func=init(), frames=np.arange(1,80), interval=1000, repeat=False)
writer = animation.FFMpegWriter(fps=1, bitrate=5000)
ani.save(workpath+"/animation_maps_GENN.mp4", writer = writer)
plt.close()
var=['Grad_gt','Grad_obs',\
     'Grad_FP_GENN_1','Grad_FP_GENN_2',\
     'Grad_FP_GENN_3','Grad_FP_GENN_4']
title=[r"$\nabla_{GT}$",'Obs (nadir+swot)',\
       r"$\nabla_{GENN}$ (all data)",r"$\nabla_{GENN}$ (masked data)",\
       r"$\nabla_{GENN}$ (all data + OI)",r"$\nabla_{GENN}$ (masked data + OI)"]
vmax=np.round(np.nanmax([np.abs(Grad_gt),np.abs(Grad_gt)]),2)
vmin=0.
cmap="viridis"
fig, ax = plt.subplots(2,3,figsize=(15,10),\
          subplot_kw=dict(projection=ccrs.PlateCarree(central_longitude=0.0)))
plt.subplots_adjust(hspace=0.5)
ani = animation.FuncAnimation(fig, animate, init_func=init(), frames=np.arange(1,80), interval=1000, repeat=False)
writer = animation.FFMpegWriter(fps=1, bitrate=5000)
ani.save(workpath+"/animation_grads_GENN.mp4", writer = writer)
plt.close()

