#!/usr/bin/env python

""" display_results.py: script to display maps, taylors diagrams and radially averaged power spectrums from MS-PB-AnDA and MS-VE-DINEOF with three different datasets (nadir / swot / nadirswot) """


__author__ = "Maxime Beauchamp"
__version__ = "0.1"
__date__ = "2019-12-10"
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

opt      = sys.argv[1]
type_obs = sys.argv[2]
domain   = sys.argv[3]
workpath = "/gpfsscratch/rech/yrf/uba22to/DINAE/"+domain+"/scores_final_"+opt+"_"+type_obs
scratchpath = '/gpfsscratch/rech/yrf/uba22to/DINAE/'+domain
if not os.path.exists(workpath):
    mk_dir_recursive(workpath)
else:
    shutil.rmtree(workpath)
    mk_dir_recursive(workpath)    

# Reload AnDA results (for GT and OI)
file_results_AnDA=scratchpath+"/resAnDA_"+opt+'_nadlag_0_'+type_obs+'/saved_path.pickle'
with open(file_results_AnDA, 'rb') as handle:
    AnDA_ssh_1, itrp_dineof = pickle.load(handle)


# Reload GE-NN results
file_results_GENN=scratchpath+"/resIA_"+opt+"_"+'nadlag_5_'+type_obs+'/FP_GENN_wwmissing_wOI/saved_path_019_FP_GENN_wwmissing.pickle'
with open(file_results_GENN, 'rb') as handle:
    itrp_FP_GENN = pickle.load(handle)[2]

			#*****************#
			# Display results #
			#*****************#

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

## Init variables for temporal analysis
nrmse_OI=np.zeros(len(AnDA_ssh_1.GT))
nrmse_Post_AnDA=np.zeros(len(AnDA_ssh_1.GT))
nrmse_VE_DINEOF=np.zeros(len(AnDA_ssh_1.GT))
nrmse_FP_GENN=np.zeros(len(AnDA_ssh_1.GT))

## Observations spatial coverage
nadswotcov=np.zeros(len(AnDA_ssh_1.GT))

# list of dates
lday1=[ datetime.strftime(datetime.strptime("2012-10-01",'%Y-%m-%d')\
                          + timedelta(days=60+i),"%Y-%m-%d") for i in range(20) ]
lday2=[ datetime.strftime(datetime.strptime("2012-10-01",'%Y-%m-%d')\
                          + timedelta(days=140+i),"%Y-%m-%d") for i in range(20) ]
lday3=[ datetime.strftime(datetime.strptime("2012-10-01",'%Y-%m-%d')\
                          + timedelta(days=220+i),"%Y-%m-%d") for i in range(20) ]
lday4=[ datetime.strftime(datetime.strptime("2012-10-01",'%Y-%m-%d')\
                          + timedelta(days=300+i),"%Y-%m-%d") for i in range(20) ]
lday = np.concatenate([lday1,lday2,lday3,lday4])
lday2 = [ datetime.strptime(lday[i],'%Y-%m-%d') for i in range(len(lday)) ] 

## Spatial analysis

for i in range(0,len(AnDA_ssh_1.GT)):

    day=lday[i]
    print(day)
    # Load data
    gt 				= AnDA_ssh_1.GT[i,:indLat,:indLon]
    Grad_gt             	= Gradient(gt,2)
    OI                          = AnDA_ssh_1.itrp_OI[i,:indLat,:indLon]
    Grad_OI                     = Gradient(OI,2)
    obs                         = AnDA_ssh_1.Obs[i,:indLat,:indLon]
    Grad_obs                    = Gradient(obs,2)
    VE_DINEOF                   = itrp_dineof[i,:indLat,:indLon]
    Grad_VE_DINEOF              = Gradient(VE_DINEOF,2)
    Post_AnDA                   = AnDA_ssh_1.itrp_postAnDA[i,:indLat,:indLon]
    Grad_Post_AnDA              = Gradient(Post_AnDA,2)
    FP_GENN                     = itrp_FP_GENN[i,:indLat,:indLon]
    Grad_FP_GENN                = Gradient(FP_GENN,2)

    ## Compute spatial coverage
    nadswotcov[i]	= len(np.argwhere(np.isfinite(obs.flatten())))/len(obs.flatten())

    ## Compute NRMSE statistics (i.e. RMSE/stdev(gt) )
    nrmse_OI[i]                  = (np.sqrt(np.nanmean(((gt-np.nanmean(gt))-(OI-np.nanmean(OI)))**2)))/np.nanstd(gt)
    nrmse_Post_AnDA[i] = (np.sqrt(np.nanmean(((gt-np.nanmean(gt))-(Post_AnDA-np.nanmean(Post_AnDA)))**2)))/np.nanstd(gt)
    nrmse_VE_DINEOF[i] = (np.sqrt(np.nanmean(((gt-np.nanmean(gt))-(VE_DINEOF-np.nanmean(VE_DINEOF)))**2)))/np.nanstd(gt)
    nrmse_FP_GENN[i]             = (np.sqrt(np.nanmean(((gt-np.nanmean(gt))-(FP_GENN-np.nanmean(FP_GENN)))**2)))/np.nanstd(gt)

    ## Display maps
    var=['obs','OI','Post_AnDA','FP_GENN','VE_DINEOF']
    title=['Obs','OI','Post-AnDA','FP-GENN','VE-DINEOF']
    fig, ax = plt.subplots(2,3,figsize=(15,10),
                          subplot_kw=dict(projection=ccrs.PlateCarree(central_longitude=0.0)))
    # display GT (reference)
    vmin = np.nanmin(gt) ; vmax = np.nanmax(gt)
    cmap="coolwarm"
    plot(ax,0,0,lon,lat,gt,'GT',\
             extent=extent,cmap=cmap,vmin=vmin,vmax=vmax)
    for ivar in range(1,len(var)+1):
        i = int(np.floor(ivar/3))
        j = int(np.floor(ivar%3))
        plot(ax,i,j,lon,lat,eval(var[ivar-1]),title[ivar-1],\
             extent=extent,cmap=cmap,vmin=vmin,vmax=vmax)
    plt.subplots_adjust(hspace=0.3,wspace=0.5)
    resfile=workpath+"/results_final_maps_"+day+".png"
    plt.savefig(resfile)       # save the figure
    plt.close()                 # close the figure

    ## Display grads
    var=['Grad_obs', 'Grad_OI','Grad_Post_AnDA','Grad_FP_GENN','Grad_VE_DINEOF']
    title=[r"$\nabla_{Obs}$",r"$\nabla_{OI}$",r"$\nabla_{Post-AnDA}$",r"$\nabla_{FP-GENN}$",r"$\nabla_{VE-DINEOF}$"]
    fig, ax = plt.subplots(2,3,figsize=(15,10),
                          subplot_kw=dict(projection=ccrs.PlateCarree(central_longitude=0.0)))
    # display GT (reference)
    vmin = np.nanmin(Grad_gt) ; vmax = np.nanmax(Grad_gt)
    cmap="viridis"
    plot(ax,0,0,lon,lat,Grad_gt,r"$\nabla_{GT}$",\
             extent=extent,cmap=cmap,vmin=vmin,vmax=vmax)
    for ivar in range(1,len(var)+1):
        i = int(np.floor(ivar/3))
        j = int(np.floor(ivar%3))
        plot(ax,i,j,lon,lat,eval(var[ivar-1]),title[ivar-1],\
             extent=extent,cmap=cmap,vmin=vmin,vmax=vmax)
    plt.subplots_adjust(hspace=0.3,wspace=0.5)
    resfile=workpath+"/results_final_grads_"+day+".png"
    plt.savefig(resfile)       # save the figure
    plt.close()                 # close the figure

## Scores tables
index=list(range(5,16))
index.extend(range(25,36))
index.extend(range(45,56))
index.extend(range(65,76))
tab_scores = np.zeros((4,3))
tab_scores[0,0] = np.nanmean(nrmse_OI[index])
tab_scores[0,1] = np.percentile(nrmse_OI[index],5)
tab_scores[0,2] = np.percentile(nrmse_OI[index],95)
tab_scores[1,0] = np.nanmean(nrmse_VE_DINEOF[index])
tab_scores[1,1] = np.percentile(nrmse_VE_DINEOF[index],5)
tab_scores[1,2] = np.percentile(nrmse_VE_DINEOF[index],95)
tab_scores[2,0] = np.nanmean(nrmse_Post_AnDA[index])
tab_scores[2,1] = np.percentile(nrmse_Post_AnDA[index],5)
tab_scores[2,2] = np.percentile(nrmse_Post_AnDA[index],95)
tab_scores[3,0] = np.nanmean(nrmse_FP_GENN[index])
tab_scores[3,1] = np.percentile(nrmse_FP_GENN[index],5)
tab_scores[3,2] = np.percentile(nrmse_FP_GENN[index],95)
np.savetxt(fname=workpath+"/tab_scores.txt",X=tab_scores,fmt='%2.2f')

## Taylor diagrams ##
# apply HPF to filter out large scale components
HR = AnDA_ssh_1.itrp_OI[:,:indLat,:indLon]
lr = np.copy(HR).reshape(HR.shape[0],-1)
tmp = lr[0,:]
sea_v2 = np.where(~np.isnan(tmp))[0]
lr_no_land = lr[:,sea_v2]
Neof=1
pca = PCA(n_components=.75)
score_global = pca.fit_transform(lr_no_land)
coeff_global = pca.components_.T
mu_global = pca.mean_
DataReconstructed_global = np.dot(score_global, coeff_global.T) +mu_global
lr[:,sea_v2] = DataReconstructed_global
lr = lr.reshape(HR.shape).flatten()
# plot Taylor diagram
resfile=workpath+"/Taylor_diagram.png"
label=['GT',\
       'OI',\
       'Post-AnDA',\
       'FP-GENN',\
       'VE-DINEOF']
series={'gt':AnDA_ssh_1.GT[:,:indLat,:indLon].flatten()-lr,
        'OI':AnDA_ssh_1.itrp_OI[:,:indLat,:indLon].flatten()-lr,
        'Post_AnDA':AnDA_ssh_1.itrp_postAnDA[:,:indLat,:indLon].flatten()-lr,
        'FP_GENN':itrp_FP_GENN[:,:indLat,:indLon].flatten()-lr,
        'VE_DINEOF':itrp_dineof[:,:indLat,:indLon].flatten()-lr}
Taylor_diag(series,label,\
            styles=['k','s','p','h','o'],\
            colors=['k','red','seagreen','darkorange','steelblue'])
plt.savefig(resfile)
plt.close()

resssh=4
## Plot averaged RAPS ##
resfile=workpath+"/results_AnDA_avg_RAPS.png"
f_ref, Pf_GT     = avg_raPsd2dv1(AnDA_ssh_1.GT[:,:indLat,:indLon],resssh,True)
f0, Pf_OI        = avg_raPsd2dv1(AnDA_ssh_1.itrp_OI[:,:indLat,:indLon],resssh,True)
f1, Pf_postAnDA  = avg_raPsd2dv1(AnDA_ssh_1.itrp_postAnDA[:,:indLat,:indLon],resssh,True)
f2, Pf_VE_DINEOF = avg_raPsd2dv1(itrp_dineof[:,:indLat,:indLon],resssh,True)
f3, Pf_FP_GENN   = avg_raPsd2dv1(itrp_FP_GENN[:,:indLat,:indLon],resssh,True)
wf_ref = 1/f_ref
wf0    = 1/f0
wf1    = 1/f1
wf2    = 1/f2
wf3    = 1/f3
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(wf_ref[1:],Pf_GT[1:],label='GT')
ax.plot(wf0[1:],Pf_OI[1:],'r-',linewidth=1.5,label='OI')
ax.plot(wf1[1:],Pf_postAnDA[1:],linestyle='solid',color='seagreen',linewidth=1.5,label='post-AnDA')
ax.plot(wf2[1:],Pf_VE_DINEOF[1:],linestyle='solid',color='steelblue',linewidth=1.5,label='VE-DINEOF')
ax.plot(wf3[1:],Pf_FP_GENN[1:],linestyle='solid',color='darkorange',linewidth=1.5,label='FP-GENN')
ax.set_xlabel("Wavenumber", fontweight='bold')
ax.set_ylabel("Power spectral density (m2/(cy/km))", fontweight='bold')
ax.set_xscale('log') ; ax.set_yscale('log')
plt.legend(loc='best',prop=dict(size='small'),frameon=False)
plt.xticks([50, 100, 200, 500, 1000], ["50km", "100km", "200km", "500km", "1000km"])
ax.invert_xaxis()
plt.grid(which='both', linestyle='--')
plt.savefig(resfile)# save the figure
plt.close() # close the figure

## Plot averaged normalize error RAPS ##
resfile=workpath+"/results_diffAnDA_avg_RAPS.png"
fig = plt.figure()
ax = fig.add_subplot(111)
f0, Pf_OI        = avg_err_raPsd2dv1(AnDA_ssh_1.itrp_OI[:,:indLat,:indLon],AnDA_ssh_1.GT[:,:indLat,:indLon],resssh,True)
f1, Pf_postAnDA  = avg_err_raPsd2dv1(AnDA_ssh_1.itrp_postAnDA[:,:indLat,:indLon],AnDA_ssh_1.GT[:,:indLat,:indLon],resssh,True)
f2, Pf_VE_DINEOF = avg_err_raPsd2dv1(itrp_dineof[:,:indLat,:indLon],AnDA_ssh_1.GT[:,:indLat,:indLon],resssh,True)
f3, Pf_FP_GENN   = avg_err_raPsd2dv1(itrp_FP_GENN[:,:indLat,:indLon],AnDA_ssh_1.GT[:,:indLat,:indLon],resssh,True)
wf0    = 1/f0
wf1    = 1/f1
wf2    = 1/f2
wf3    = 1/f3
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(wf0[1:],Pf_OI[1:],'r-',linewidth=1.5,label='OI')
ax.plot(wf1[1:],Pf_postAnDA[1:],linestyle='solid',color='seagreen',linewidth=1.5,label='post-AnDA')
ax.plot(wf2[1:],Pf_VE_DINEOF[1:],linestyle='solid',color='steelblue',linewidth=1.5,label='VE-DINEOF')
ax.plot(wf3[1:],Pf_FP_GENN[1:],linestyle='solid',color='darkorange',linewidth=1.5,label='FP-GENN')
ax.set_xlabel("Wavenumber", fontweight='bold')
ax.set_ylabel("Signal-to-noise ratio", fontweight='bold')
ax.set_xscale('log') ; ax.set_yscale('log')
plt.legend(loc='best',prop=dict(size='small'),frameon=False)
plt.xticks([50, 100, 200, 500, 1000], ["50km", "100km", "200km", "500km", "1000km"])
ax.invert_xaxis()
plt.grid(which='both', linestyle='--')
plt.savefig(resfile)# save the figure
plt.close() # close the figure

## NRMSE (time series) ##
N = len(lday)
print(N)
# first axis with nRMSE time series
plt.plot(range(N),nrmse_OI,linestyle='solid',color='red',linewidth=1.5,label='OI')
plt.plot(range(N),nrmse_Post_AnDA,linestyle='solid',color='seagreen',linewidth=1.5,label='post-AnDA')
plt.plot(range(N),nrmse_VE_DINEOF,linestyle='solid',color='steelblue',linewidth=1.5,label='VE-DINEOF')
plt.plot(range(N),nrmse_FP_GENN,linestyle='solid',color='darkorange',linewidth=1.5,label='FP-GENN')
# add vertical bar to divide the 4 periods
plt.axvline(x=19)
plt.axvline(x=39)
plt.axvline(x=59)
# graphical options
plt.ylabel('nRMSE')
plt.xlabel('Time (days)')
plt.xticks([0,20,40,60],\
           [lday[0],lday[20],lday[40],lday[60]],\
           rotation=45, ha='right')
plt.margins(x=0)
plt.grid(True,alpha=.3)
plt.legend(loc='upper left',prop=dict(size='small'),frameon=False,bbox_to_anchor=(0,1.02,1,0.2),ncol=2,mode="expand")
# second axis with spatial coverage
axes2 = plt.twinx()
width=0.75
p1 = axes2.bar(range(N), nadswotcov, width,color='r',alpha=0.25)
axes2.set_ylim(0, 1)
axes2.set_ylabel('Spatial Coverage (%)')
axes2.margins(x=0)
resfile=workpath+"/TS_AnDA_nRMSE.png"
plt.savefig(resfile,bbox_inches="tight")    # save the figure
plt.close()         	# close the figure

## animations (videos) ##
def init():
    global fig, ax
    # Load data
    gt                          = AnDA_ssh_1.GT[0,:indLat,:indLon]
    Grad_gt                     = Gradient(gt,2)
    OI                          = AnDA_ssh_1.itrp_OI[0,:indLat,:indLon]
    Grad_OI                     = Gradient(OI,2)
    obs                         = AnDA_ssh_1.Obs[0,:indLat,:indLon]
    Grad_obs                    = Gradient(obs,2)
    Post_AnDA                   = AnDA_ssh_1.itrp_postAnDA[0,:indLat,:indLon]
    Grad_Post_AnDA              = Gradient(Post_AnDA,2)
    FP_GENN                     = itrp_FP_GENN[0,:indLat,:indLon]
    Grad_FP_GENN                = Gradient(FP_GENN,2)
    VE_DINEOF                   = itrp_dineof[0,:indLat,:indLon]
    Grad_VE_DINEOF              = Gradient(VE_DINEOF,2)
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
    OI                          = AnDA_ssh_1.itrp_OI[i,:indLat,:indLon]
    Grad_OI                     = Gradient(OI,2)
    obs                         = AnDA_ssh_1.Obs[i,:indLat,:indLon]
    Grad_obs                    = Gradient(obs,2)
    Post_AnDA                   = AnDA_ssh_1.itrp_postAnDA[i,:indLat,:indLon]
    Grad_Post_AnDA              = Gradient(Post_AnDA,2)
    FP_GENN                     = itrp_FP_GENN[i,:indLat,:indLon]
    Grad_FP_GENN                = Gradient(FP_GENN,2)
    VE_DINEOF                   = itrp_dineof[i,:indLat,:indLon]
    Grad_VE_DINEOF              = Gradient(VE_DINEOF,2)
    for ivar in range(len(var)):
        ii = int(np.floor(ivar/3))
        jj = int(np.floor(ivar%3))
        ax[ii][jj].cla()
        plot(ax,ii,jj,lon,lat,eval(var[ivar]),title[ivar],\
              extent=extent,cmap=cmap,vmin=vmin,vmax=vmax,colorbar=False)
    return fig, ax

var=['gt','obs',\
     'OI','Post_AnDA',\
     'FP_GENN','VE_DINEOF']
title=['Ground Truth','Obs',\
       'OI','Post-AnDA',\
       'FP_GENN','VE-DINEOF']
vmax=np.round(np.nanmax([np.abs(AnDA_ssh_1.GT),np.abs(AnDA_ssh_1.GT)]),2)
vmin=-1.0*vmax
cmap="coolwarm"
fig, ax = plt.subplots(2,3,figsize=(15,10),\
          subplot_kw=dict(projection=ccrs.PlateCarree(central_longitude=0.0)))
plt.subplots_adjust(hspace=0.5)
ani = animation.FuncAnimation(fig, animate, init_func=init(), frames=np.arange(1,80), interval=1000, repeat=False)
writer = animation.FFMpegWriter(fps=1, bitrate=5000)
ani.save(workpath+"/animation_maps.mp4", writer = writer)
plt.close()
var=['Grad_gt','Grad_obs',\
     'Grad_OI','Grad_Post_AnDA',\
     'Grad_FP_GENN','Grad_VE_DINEOF']
title=[r"$\nabla_{GT}$",'Obs',\
       r"$\nabla_{OI}$",r"$\nabla_{Post-AnDA}$",\
       r"$\nabla_{FP-GENN}$",r"$\nabla_{VE-DINEOF}$"]
vmax=np.round(np.nanmax([np.abs(Grad_gt),np.abs(Grad_gt)]),2)
vmin=0.
cmap="viridis"
fig, ax = plt.subplots(2,3,figsize=(15,10),\
          subplot_kw=dict(projection=ccrs.PlateCarree(central_longitude=0.0)))
plt.subplots_adjust(hspace=0.5)
ani = animation.FuncAnimation(fig, animate, init_func=init(), frames=np.arange(1,80), interval=1000, repeat=False)
writer = animation.FFMpegWriter(fps=1, bitrate=5000)
ani.save(workpath+"/animation_grads.mp4", writer = writer)
plt.close()

