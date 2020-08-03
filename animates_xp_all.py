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

AnDA_lag = sys.argv[1]
NN_lag   = sys.argv[2]
type_obs = sys.argv[3]
domain   = sys.argv[4] 
workpath = "/gpfsscratch/rech/yrf/uba22to/DINAE/"+domain+"/animates_allmethods_AnDAnadlag_"+AnDA_lag+"_NNnadlag_"+NN_lag+"_"+type_obs
scratchpath = '/gpfsscratch/rech/yrf/uba22to/DINAE/'+domain
if not os.path.exists(workpath):
    mk_dir_recursive(workpath)
else:
    shutil.rmtree(workpath)
    mk_dir_recursive(workpath)

# Reload saved AnDA result
file_results_nadir=scratchpath+'/resAnDA_nadir_nadlag_'+AnDA_lag+"_"+type_obs+'/saved_path.pickle'
with open(file_results_nadir, 'rb') as handle:
    AnDA_ssh_1, itrp_dineof = pickle.load(handle)
    AnDA_ssh_1_nadir = AnDA_ssh_1
    itrp_dineof_nadir = itrp_dineof
'''file_results_swot=scratchpath+'/resAnDA_swot_nadlag_0_'+type_obs+'/saved_path.pickle'
with open(file_results_swot, 'rb') as handle:
    AnDA_ssh_1, itrp_dineof = pickle.load(handle)
    AnDA_ssh_1_swot = AnDA_ssh_1
    itrp_dineof_swot = itrp_dineof'''
file_results_nadirswot=scratchpath+'/resAnDA_nadirswot_nadlag_'+AnDA_lag+"_"+type_obs+'/saved_path.pickle'
with open(file_results_nadirswot, 'rb') as handle:
    AnDA_ssh_1, itrp_dineof = pickle.load(handle)
    AnDA_ssh_1_nadirswot = AnDA_ssh_1
    itrp_dineof_nadirswot = itrp_dineof

# Reload saved ConvAE and GE-NN results
file_results_nadir=scratchpath+'/resIA_nadir_nadlag_'+NN_lag+"_"+type_obs+'/FP_ConvAE_womissing_wocov/saved_path_019_FP_ConvAE_womissing.pickle'
with open(file_results_nadir, 'rb') as handle:
    itrp_FP_ConvAE_nadir = pickle.load(handle)[2]
file_results_nadirswot=scratchpath+'/resIA_nadirswot_nadlag_'+NN_lag+"_"+type_obs+'/FP_ConvAE_womissing_wocov/saved_path_019_FP_ConvAE_womissing.pickle'
with open(file_results_nadirswot, 'rb') as handle:
    itrp_FP_ConvAE_nadirswot = pickle.load(handle)[2]
file_results_nadir=scratchpath+'/resIA_nadir_nadlag_'+NN_lag+"_"+type_obs+'/FP_GENN_wmissing_wOI/saved_path_019_FP_GENN_wmissing.pickle'
with open(file_results_nadir, 'rb') as handle:
    itrp_FP_GENN_nadir = pickle.load(handle)[2]
file_results_nadirswot=scratchpath+'/resIA_nadirswot_nadlag_'+NN_lag+"_"+type_obs+'/FP_GENN_wmissing_wOI/saved_path_019_FP_GENN_wmissing.pickle'
with open(file_results_nadirswot, 'rb') as handle:
    itrp_FP_GENN_nadirswot = pickle.load(handle)[2]

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
Grad_gt=np.zeros((len(AnDA_ssh_1.GT),indLat,indLon))
for i in range(len(AnDA_ssh_1.GT)):
    Grad_gt[i,:,:] = Gradient(AnDA_ssh_1.GT[i,:indLat,:indLon],2)

def init(tobs):
    print(tobs)
    global fig, ax
    # Load data
    gt                          = AnDA_ssh_1.GT[0,:indLat,:indLon]
    Grad_gt                     = Gradient(gt,2)
    OI                          = AnDA_ssh_1.itrp_OI[0,:indLat,:indLon]
    Grad_OI                     = Gradient(OI,2)
    obs                         = eval("AnDA_ssh_1_"+tobs).Obs[0,:indLat,:indLon]
    Grad_obs                    = Gradient(obs,2)
    Post_AnDA                   = eval("AnDA_ssh_1_"+tobs).itrp_postAnDA[0,:indLat,:indLon]
    Grad_Post_AnDA              = Gradient(Post_AnDA,2)
    FP_ConvAE                   = eval("itrp_FP_ConvAE_"+tobs)[0,:indLat,:indLon]
    Grad_FP_ConvAE              = Gradient(FP_ConvAE,2)
    FP_GENN                     = eval("itrp_FP_GENN_"+tobs)[0,:indLat,:indLon]
    Grad_FP_GENN                = Gradient(FP_GENN,2)
    for ivar in range(len(var)):
        ii = int(np.floor(ivar/3))
        jj = int(np.floor(ivar%3))
        plot(ax,ii,jj,lon,lat,eval(var[ivar]),title[ivar],\
              extent=extent,cmap=cmap,vmin=vmin,vmax=vmax,colorbar=True,\
              orientation="vertical")

def animate(i,tobs):
    global fig, ax
    # Load data
    gt                          = AnDA_ssh_1.GT[i,:indLat,:indLon]
    Grad_gt                     = Gradient(gt,2)
    OI                          = AnDA_ssh_1.itrp_OI[i,:indLat,:indLon]
    Grad_OI                     = Gradient(OI,2)
    obs                         = eval("AnDA_ssh_1_"+tobs).Obs[i,:indLat,:indLon]
    Grad_obs                    = Gradient(obs,2)
    Post_AnDA                   = eval("AnDA_ssh_1_"+tobs).itrp_postAnDA[i,:indLat,:indLon]
    Grad_Post_AnDA              = Gradient(Post_AnDA,2)
    FP_ConvAE                   = eval("itrp_FP_ConvAE_"+tobs)[i,:indLat,:indLon]
    Grad_FP_ConvAE              = Gradient(FP_ConvAE,2)
    FP_GENN                     = eval("itrp_FP_GENN_"+tobs)[i,:indLat,:indLon]
    Grad_FP_GENN                = Gradient(FP_GENN,2)
    for ivar in range(len(var)):   
        ii = int(np.floor(ivar/3))
        jj = int(np.floor(ivar%3))
        ax[ii][jj].cla()
        plot(ax,ii,jj,lon,lat,eval(var[ivar]),title[ivar],\
              extent=extent,cmap=cmap,vmin=vmin,vmax=vmax,colorbar=False)
    return fig, ax

# for nadir
var=['gt','obs',\
     'OI','Post_AnDA',\
     'FP_ConvAE','FP_GENN']
title=['Ground Truth','Obs (nadir)',\
       'OI (nadir)','Post-AnDA (nadir)',\
       'FP_ConvAE (nadir)',\
       'FP_GENN (nadir)']
vmax=np.round(np.nanmax([np.abs(AnDA_ssh_1.GT),np.abs(AnDA_ssh_1.GT)]),2)
vmin=-1.0*vmax
cmap="coolwarm"
fig, ax = plt.subplots(2,3,figsize=(15,10),\
          subplot_kw=dict(projection=ccrs.PlateCarree(central_longitude=0.0)))
plt.subplots_adjust(hspace=0.5)
ani = animation.FuncAnimation(fig, animate, init_func=init("nadir"), frames=np.arange(1,80), fargs=("nadir",), interval=1000, repeat=False)
writer = animation.FFMpegWriter(fps=1, bitrate=5000)
ani.save(workpath+"/animation_maps_nadir.mp4", writer = writer)
plt.close()
var=['Grad_gt','Grad_obs',\
     'Grad_OI','Grad_Post_AnDA',\
     'Grad_FP_ConvAE','Grad_FP_GENN']
title=[r"$\nabla_{GT}$",'Obs (nadir)',\
       r"$\nabla_{OI}$ (nadir)",r"$\nabla_{Post-AnDA}$ (nadir)",\
       r"$\nabla_{FP-ConvAE}$ (nadir)",r"$\nabla_{FP-GENN}$ (nadir)"]
vmax=np.round(np.nanmax([np.abs(Grad_gt),np.abs(Grad_gt)]),2)
vmin=0.
cmap="viridis"
fig, ax = plt.subplots(2,3,figsize=(15,10),\
          subplot_kw=dict(projection=ccrs.PlateCarree(central_longitude=0.0)))
plt.subplots_adjust(hspace=0.5)
ani = animation.FuncAnimation(fig, animate, init_func=init("nadir"), frames=np.arange(1,80), fargs=("nadir",), interval=1000, repeat=False)
writer = animation.FFMpegWriter(fps=1, bitrate=5000)
ani.save(workpath+"/animation_grads_nadir.mp4", writer = writer)
plt.close()

# for nadirswot
var=['gt','obs',\
     'OI','Post_AnDA',\
     'FP_ConvAE','FP_GENN']
title=['Ground Truth','Obs (nadir+swot)',\
       'OI (nadir)','Post-AnDA (nadir+swot)',\
       'FP_ConvAE (nadir+swot)',\
       'FP_GENN (nadir+swot)']
vmax=np.round(np.nanmax([np.abs(AnDA_ssh_1.GT),np.abs(AnDA_ssh_1.GT)]),2)
vmin=-1.0*vmax
cmap="coolwarm"
fig, ax = plt.subplots(2,3,figsize=(15,10),\
          subplot_kw=dict(projection=ccrs.PlateCarree(central_longitude=0.0))) 
plt.subplots_adjust(hspace=0.5)
ani = animation.FuncAnimation(fig, animate, init_func=init("nadirswot"), frames=np.arange(1,80), fargs=("nadirswot",), interval=1000, repeat=False)
writer = animation.FFMpegWriter(fps=1, bitrate=5000)
ani.save(workpath+"/animation_maps_nadirswot.mp4", writer = writer)
plt.close()
var=['Grad_gt','Grad_obs',\
     'Grad_OI','Grad_Post_AnDA',\
     'Grad_FP_ConvAE','Grad_FP_GENN']
title=[r"$\nabla_{GT}$",'Obs (nadir+swot)',\
       r"$\nabla_{OI}$ (nadir)",r"$\nabla_{Post-AnDA}$ (nadir+swot)",\
       r"$\nabla_{FP-ConvAE}$ (nadir+swot)",r"$\nabla_{FP-GENN}$ (nadir+swot)"]
vmax=np.round(np.nanmax([np.abs(Grad_gt),np.abs(Grad_gt)]),2)
vmin=0.
cmap="viridis"
fig, ax = plt.subplots(2,3,figsize=(15,10),\
          subplot_kw=dict(projection=ccrs.PlateCarree(central_longitude=0.0)))
plt.subplots_adjust(hspace=0.5)
ani = animation.FuncAnimation(fig, animate, init_func=init("nadirswot"), frames=np.arange(1,80), fargs=("nadirswot",), interval=1000, repeat=False)
writer = animation.FFMpegWriter(fps=1, bitrate=5000)
ani.save(workpath+"/animation_grads_nadirswot.mp4", writer = writer)
plt.close()

