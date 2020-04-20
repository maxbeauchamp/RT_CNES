#!/usr/bin/env python

""" display_results.py: script to display maps, taylors diagrams and radially averaged power spectrums from MS-PB-AnDA and MS-VE-DINEOF with three different datasets (nadir / swot / nadirswot) """

__author__ = "Maxime Beauchamp"
__version__ = "0.1"
__date__ = "2019-12-10"
__email__ = "maxime.beauchamp76@gmail.com"

from pb_anda import *
import matplotlib.dates as mdates
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
workpath = "/home3/scratch/mbeaucha/scores_allmethods_AnDAnadlag_"+AnDA_lag+"_NNnadlag_"+NN_lag+"_"+type_obs
if not os.path.exists(workpath):
    mk_dir_recursive(workpath)
else:
    shutil.rmtree(workpath)
    mk_dir_recursive(workpath)    

# Reload saved AnDA result
file_results_nadir='/home3/scratch/mbeaucha/resAnDA_nadir_nadlag_'+AnDA_lag+"_"+type_obs+'/saved_path.pickle'
with open(file_results_nadir, 'rb') as handle:
    AnDA_ssh_1, itrp_dineof = pickle.load(handle)
    AnDA_ssh_1_nadir = AnDA_ssh_1  
    itrp_dineof_nadir = itrp_dineof
file_results_nadirswot='/home3/scratch/mbeaucha/resAnDA_nadirswot_nadlag_'+AnDA_lag+"_"+type_obs+'/saved_path.pickle'
with open(file_results_nadirswot, 'rb') as handle:
    AnDA_ssh_1, itrp_dineof = pickle.load(handle)
    AnDA_ssh_1_nadirswot = AnDA_ssh_1  
    itrp_dineof_nadirswot = itrp_dineof

# Reload saved ConvAE and GE-NN results
file_results_nadir='/home3/scratch/mbeaucha/resIA_nadir_nadlag_'+NN_lag+"_"+type_obs+'/FP_ConvAE_womissing_wocov/saved_path_019_FP_ConvAE_womissing.pickle'
with open(file_results_nadir, 'rb') as handle:
    itrp_FP_ConvAE_nadir = pickle.load(handle)[2]
file_results_nadirswot='/home3/scratch/mbeaucha/resIA_nadirswot_nadlag_'+NN_lag+"_"+type_obs+'/FP_ConvAE_womissing_wocov/saved_path_019_FP_ConvAE_womissing.pickle'
with open(file_results_nadirswot, 'rb') as handle:
    itrp_FP_ConvAE_nadirswot = pickle.load(handle)[2]
file_results_nadir='/home3/scratch/mbeaucha/resIA_nadir_nadlag_'+NN_lag+"_"+type_obs+'/FP_GENN_womissing_wocov/saved_path_019_FP_GENN_womissing.pickle'
with open(file_results_nadir, 'rb') as handle:
    itrp_FP_GENN_nadir = pickle.load(handle)[2]
file_results_nadirswot='/home3/scratch/mbeaucha/resIA_nadirswot_nadlag_'+NN_lag+"_"+type_obs+'/FP_GENN_womissing_wocov/saved_path_019_FP_GENN_womissing.pickle'
with open(file_results_nadirswot, 'rb') as handle:
    itrp_FP_GENN_nadirswot = pickle.load(handle)[2]

			#*****************#
			# Display results #
			#*****************#
resssh = 4
lon = np.arange(-65,-55,1/20)
lat = np.arange(30,40,1/20)
indLat  = np.arange(0,200)
indLon  = np.arange(0,200)
lon = lon[indLon]
lat = lat[indLat]
extent_ = [np.min(lon),np.max(lon),np.min(lat),np.max(lat)]

## Init variables for temporal analysis
nrmse_OI=np.zeros(len(AnDA_ssh_1.GT))
nrmse_Post_AnDA_nadir=np.zeros(len(AnDA_ssh_1.GT))
nrmse_VE_DINEOF_nadir=np.zeros(len(AnDA_ssh_1.GT))
nrmse_FP_ConvAE_nadir=np.zeros(len(AnDA_ssh_1.GT))
nrmse_FP_GENN_nadir=np.zeros(len(AnDA_ssh_1.GT))
nrmse_Post_AnDA_nadirswot=np.zeros(len(AnDA_ssh_1.GT))
nrmse_VE_DINEOF_nadirswot=np.zeros(len(AnDA_ssh_1.GT))
nrmse_FP_ConvAE_nadirswot=np.zeros(len(AnDA_ssh_1.GT))
nrmse_FP_GENN_nadirswot=np.zeros(len(AnDA_ssh_1.GT))
## Observations spatial coverage
nadcov=np.zeros(len(AnDA_ssh_1.GT))
swotcov=np.zeros(len(AnDA_ssh_1.GT))
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

indLon=200
indLat=200

for i in range(0,len(AnDA_ssh_1.GT)):

    day=lday[i]
    print(day)
    ## Maps
    resfile1=workpath+"/results_AnDA_maps_"+day+".png"
    resfile2=workpath+"/results_AnDA_grads_"+day+".png"
    # Load data
    gt 				= AnDA_ssh_1.GT[i,:indLon,:indLat]
    Grad_gt             	= Gradient(gt,2)
    # nadir
    OI_nadir                    = AnDA_ssh_1_nadir.itrp_OI[i,:indLon,:indLat]
    Grad_OI_nadir               = Gradient(OI_nadir,2)
    obs_nadir 			= AnDA_ssh_1_nadir.Obs[i,:indLon,:indLat]
    VE_DINEOF_nadir           	= itrp_dineof_nadir[i,:indLon,:indLat]
    Grad_VE_DINEOF_nadir 	= Gradient(VE_DINEOF_nadir,2)
    AnDA_nadir 		 	= AnDA_ssh_1_nadir.itrp_AnDA[i,:indLon,:indLat]
    Grad_AnDA_nadir          	= Gradient(AnDA_nadir,2)
    Post_AnDA_nadir 		= AnDA_ssh_1_nadir.itrp_postAnDA[i,:indLon,:indLat]
    Grad_Post_AnDA_nadir	= Gradient(Post_AnDA_nadir,2)
    FP_ConvAE_nadir		= itrp_FP_ConvAE_nadir[i,:indLon,:indLat]
    Grad_FP_ConvAE_nadir        = Gradient(FP_ConvAE_nadir,2)
    FP_GENN_nadir               = itrp_FP_GENN_nadir[i,:indLon,:indLat]
    Grad_FP_GENN_nadir          = Gradient(FP_GENN_nadir,2)
    # nadirswot
    OI_nadirswot                = AnDA_ssh_1_nadirswot.itrp_OI[i,:indLon,:indLat]
    Grad_OI_nadirswot           = Gradient(OI_nadirswot,2)
    obs_nadirswot 		= AnDA_ssh_1_nadirswot.Obs[i,:indLon,:indLat]
    VE_DINEOF_nadirswot         = itrp_dineof_nadirswot[i,:indLon,:indLat]
    Grad_VE_DINEOF_nadirswot    = Gradient(VE_DINEOF_nadirswot,2)
    AnDA_nadirswot 		= AnDA_ssh_1_nadirswot.itrp_AnDA[i,:indLon,:indLat]
    Grad_AnDA_nadirswot         = Gradient(AnDA_nadirswot,2)
    Post_AnDA_nadirswot         = AnDA_ssh_1_nadirswot.itrp_postAnDA[i,:indLon,:indLat]
    Grad_Post_AnDA_nadirswot    = Gradient(Post_AnDA_nadirswot,2)
    FP_ConvAE_nadirswot         = itrp_FP_ConvAE_nadirswot[i,:indLon,:indLat]
    Grad_FP_ConvAE_nadirswot    = Gradient(FP_ConvAE_nadirswot,2)
    FP_GENN_nadirswot           = itrp_FP_GENN_nadirswot[i,:indLon,:indLat]
    Grad_FP_GENN_nadirswot      = Gradient(FP_GENN_nadirswot,2)

    ## Compute spatial coverage
    nadcov[i]		= len(np.argwhere(np.isfinite(obs_nadir.flatten())))/len(obs_nadir.flatten())
    swotcov[i]		= len(np.argwhere(np.isfinite(obs_swot.flatten())))/len(obs_swot.flatten())
    nadswotcov[i]	= len(np.argwhere(np.isfinite(obs_nadirswot.flatten())))/len(obs_nadirswot.flatten())

    ## Compute NRMSE statistics (i.e. RMSE/stdev(gt) )
    nrmse_OI_nadir[i]		= (np.sqrt(np.nanmean(((gt-np.nanmean(gt))-(OI_nadir-np.nanmean(OI_nadir)))**2)))/np.nanstd(gt)
    nrmse_Post_AnDA_nadir[i]	= (np.sqrt(np.nanmean(((gt-np.nanmean(gt))-(Post_AnDA_nadir-np.nanmean(Post_AnDA_nadir)))**2)))/np.nanstd(gt)
    nrmse_VE_DINEOF_nadir[i]	= (np.sqrt(np.nanmean(((gt-np.nanmean(gt))-(VE_DINEOF_nadir-np.nanmean(VE_DINEOF_nadir)))**2)))/np.nanstd(gt)
    nrmse_FP_ConvAE_nadir[i]    = (np.sqrt(np.nanmean(((gt-np.nanmean(gt))-(FP_ConvAE_nadir-np.nanmean(FP_ConvAE_nadir)))**2)))/np.nanstd(gt)
    nrmse_FP_GENN_nadir[i]      = (np.sqrt(np.nanmean(((gt-np.nanmean(gt))-(FP_GENN_nadir-np.nanmean(FP_GENN_nadir)))**2)))/np.nanstd(gt)
    nrmse_Post_AnDA_nadirswot[i]= (np.sqrt(np.nanmean(((gt-np.nanmean(gt))-(Post_AnDA_nadirswot-np.nanmean(Post_AnDA_nadirswot)))**2)))/np.nanstd(gt)
    nrmse_VE_DINEOF_nadirswot[i]= (np.sqrt(np.nanmean(((gt-np.nanmean(gt))-(VE_DINEOF_nadirswot-np.nanmean(VE_DINEOF_nadirswot)))**2)))/np.nanstd(gt) 
    nrmse_FP_ConvAE_nadirswot[i]  = (np.sqrt(np.nanmean(((gt-np.nanmean(gt))-(FP_ConvAE_nadirswot-np.nanmean(FP_ConvAE_nadirswot)))**2)))/np.nanstd(gt)
    nrmse_FP_GENN_nadirswot[i]    = (np.sqrt(np.nanmean(((gt-np.nanmean(gt))-(FP_GENN_nadirswot-np.nanmean(FP_GENN_nadirswot)))**2)))/np.nanstd(gt)

    # Display individual maps
    var=['gt','obs_nadir','obs_nadirswot',\
         'OI_nadir','OI_nadirswot'\
         'AnDA_nadir','AnDA_nadirswot',\
         'Post_AnDA_nadir','Post_AnDA_nadirswot',\
         'FP_ConvAE_nadir','FP_ConvAE_nadirswot',\
         'FP_GENN_nadir','FP_GENN_nadirswot',\
         'VE_DINEOF_nadir','VE_DINEOF_nadirswot',]
    title=['Ground Truth','Obs (nadir)','Obs (nadir+swot)',\
           'OI (nadir)','OI (nadir+swot)',\
           'AnDA (nadir)','AnDA (nadir+swot)',\
           'Post-AnDA (nadir)','Post-AnDA (nadir+swot)',\
           'FP-ConvAE (nadir)''FP-ConvAE (nadir+swot)',\
           'FP-GENN (nadir)','FP-GENN (nadir+swot)',\
           'VE-DINEOF (nadir)','VE-DINEOF (nadir+swot)']
    for ivar in range(0,len(var)):
        resfile = workpath+"/results_AnDA_maps_"+var[ivar]+'_'+day+".png"
        fig, ax = plt.subplots(1,1,figsize=(10,10),
                          subplot_kw=dict(projection=ccrs.PlateCarree(central_longitude=0.0)))   
        vmin = -2 ; vmax = 2
        cmap="coolwarm"
        plot2(ax,lon,lat,eval(var[ivar]),title[ivar],\
             extent=extent_,cmap=cmap,vmin=vmin,vmax=vmax)
        plt.savefig(resfile)       # save the figure
        plt.close()                 # close the figure

    # Display individual gradient maps
    var=['Grad_gt',\
         'Grad_OI_nadirswot','Grad_OI_nadirswot',\
         'Grad_AnDA_nadir','Grad_AnDA_nadirswot',\
         'Grad_Post_AnDA_nadir','Grad_Post_AnDA_nadirswot',\
         'Grad_FP_ConvAE_nadir','Grad_FP_ConvAE_nadirswot',\
         'Grad_FP_GENN_nadir','Grad_FP_GENN_nadirswot',\
         'Grad_VE_DINEOF_nadir','Grad_VE_DINEOF_nadirswot']
    title=[r"$\nabla_{GT}$",\
           r"$\nabla_{OI}$ (nadir)",r"$\nabla_{OI}$ (nadir+swot)",\
           r"$\nabla_{AnDA}$ (nadir)",r"$\nabla_{AnDA} (nadir+swot)$",\
           r"$\nabla_{Post-AnDA}$ (nadir)",r"$\nabla_{Post-AnDA} (nadir+swot)$",\
           r"$\nabla_{FP-ConvAE}$ (nadir)",r"$\nabla_{FP-ConvAE} (nadir+swot)$",\
           r"$\nabla_{FP-GENN}$ (nadir)",r"$\nabla_{FP-GENN} (nadir+swot)$",\
           r"$\nabla_{VE-DINEOF}$ (nadir)",r"$\nabla_{VE-DINEOF} (nadir+swot)$",]
    for ivar in range(0,len(var)):
        resfile = workpath+"/results_AnDA_grads_"+var[ivar]+'_'+day+".png"
        fig, ax = plt.subplots(1,1,figsize=(10,10),
                          subplot_kw=dict(projection=ccrs.PlateCarree(central_longitude=0.0)))
        vmin = 0 ; vmax = 0.2
        cmap="viridis"
        plot2(ax,lon,lat,eval(var[ivar]),title[ivar],\
             extent=extent_,cmap=cmap,vmin=vmin,vmax=vmax)
        plt.savefig(resfile)       # save the figure
        plt.close()                 # close the figure

    # Display maps
    var=['obs_nadir','obs_nadirswot',\
         'OI_nadir','OI_nadirswot',\
         'AnDA_nadir','AnDA_nadirswot',\
         'Post_AnDA_nadir','Post_AnDA_nadirswot',\
         'FP_ConvAE_nadir','FP_ConvAE_nadirswot',\
         'FP_GENN_nadir','FP_GENN_nadirswot',\
         'VE_DINEOF_nadir','VE_DINEOF_nadirswot',]
    title=['Obs (nadir)','Obs (nadir+swot)',\
           'OI','OI',\
           'AnDA','AnDA',\
           'Post-AnDA','Post-AnDA',\
           'FP-ConvAE','FP-ConvAE',\
           'FP-GENN','FP-GENN',\
           'VE-DINEOF','VE-DINEOF']
    fig, ax = plt.subplots(2,8,figsize=(40,25),
                          subplot_kw=dict(projection=ccrs.PlateCarree(central_longitude=0.0)))
    # display GT (reference)
    vmin = np.nanmin(gt) ; vmax = np.nanmax(gt)
    cmap="coolwarm"
    plot(ax,0,0,lon,lat,gt,'GT',\
             extent=extent_,cmap=cmap,vmin=vmin,vmax=vmax)
    ax[1][0].set_visible(False)
    ax[2][0].set_visible(False)
    for ivar in range(0,len(var)):
        i = int(np.floor(ivar/2))+1 ; j = (ivar%2)
        plot(ax,j,i,lon,lat,eval(var[ivar]),title[ivar],\
             extent=extent_,cmap=cmap,vmin=vmin,vmax=vmax)
    plt.subplots_adjust(hspace=0.3,wspace=0.5)
    plt.savefig(resfile1)	# save the figure
    plt.close()			# close the figure

    # Display gradients
    var=['obs_nadir','obs_nadirswot',\
         'Grad_OI','Grad_OI',\
         'Grad_AnDA_nadir','Grad_AnDA_nadirswot',\
         'Grad_Post_AnDA_nadir','Grad_Post_AnDA_nadirswot',\
         'Grad_FP_ConvAE_nadir','Grad_FP_ConvAE_nadirswot',\
         'Grad_FP_GENN_nadir','Grad_FP_GENN_nadirswot',\
         'Grad_VE_DINEOF_nadir','Grad_VE_DINEOF_nadirswot']
    title=['Obs (nadir)','Obs (nadir+swot)',\
           r"$\nabla_{OI}$",r"$\nabla_{OI}$",\
           r"$\nabla_{AnDA}",r"$\nabla_{AnDA}$",\
           r"$\nabla_{Post-AnDA}$",r"$\nabla_{Post-AnDA}$",\
           r"$\nabla_{FP-ConvAE}$",r"$\nabla_{FP-ConvAE}$",\
           r"$\nabla_{FP-GENN}$",r"$\nabla_{FP-GENN}$",\
           r"$\nabla_{VE-DINEOF}$",r"$\nabla_{VE-DINEOF}$",]
    fig, ax = plt.subplots(2,8,figsize=(40,25),
                          subplot_kw=dict(projection=ccrs.PlateCarree(central_longitude=0.0)))
    vmin = np.nanmin(Grad_gt) ; vmax = np.nanmax(Grad_gt)
    vmin2 = np.nanmin(gt) ; vmax2 = np.nanmax(gt)
    cmap="viridis"
    cmap2="coolwarm"
    plot(ax,0,0,lon,lat,Grad_gt,r"$\nabla_{GT}$",\
             extent=extent_,cmap=cmap,vmin=vmin,vmax=vmax)
    ax[1][0].set_visible(False)
    ax[2][0].set_visible(False)
    for ivar in range(0,len(var)):
        i = int(np.floor(ivar/2))+1 ; j = (ivar%2)
        if i==1:
            plot(ax,j,i,lon,lat,eval(var[ivar]),title[ivar],\
                 extent=extent_,cmap=cmap2,vmin=vmin2,vmax=vmax2)
        else:
            plot(ax,j,i,lon,lat,eval(var[ivar]),title[ivar],\
                 extent=extent_,cmap=cmap,vmin=vmin,vmax=vmax)
    plt.subplots_adjust(hspace=0.3,wspace=0.5)
    plt.savefig(resfile2)	# save the figure
    plt.close()			# close the figure

    ## Taylor diagrams (rough variables)
    resfile=workpath+"/Taylor_diagram_maps_"+day+".png"
    label=['GT',\
           'OI (nadir)','OI (nadir+swot)',\
           'AnDA (nadir)','AnDA (nadir+swot)',\
           'Post-AnDA (nadir)','Post-AnDA (nadir+swot)',\
           'FP-ConvAE (nadir)','FP-ConvAE (nadir+swot)',\
           'FP-GENN (nadir)','FP-GENN (nadir+swot)',\
           'VE-DINEOF (nadir)','VE-DINEOF (nadir+swot)']
    series={'gt':gt,
            'OI_nadir':OI_nadir,
            'OI_nadirswot':OI_nadirswot,
            'AnDA_nadir':AnDA_nadir,'AnDA_nadirswot':AnDA_nadirswot,
            'Post_AnDA_nadir':Post_AnDA_nadir,'Post_AnDA_nadirswot':Post_AnDA_nadirswot,
            'FP_ConvAE_nadir':FP_ConvAE_nadir,'FP_ConvAE_nadirswot':FP_ConvAE_nadirswot,
            'FP_GENN_nadir':FP_GENN_nadir,'FP_GENN_nadirswot':FP_GENN_nadirswot,
            'VE_DINEOF_nadir':VE_DINEOF_nadir,'VE_DINEOF_nadirswot':VE_DINEOF_nadirswot}
    Taylor_diag(series,label,\
                styles=['k','s','p','o','p','o','p','o','p','o','p','o'],\
                colors=['k','y','mediumseagreen','mediumseagreen',\
                                'seagreen','seagreen',\
                                'steelblue','steelblue',\
                                'mediumorchid','mediumorchid',\
                                'darkorange','darkorange'])
    plt.savefig(resfile)
    plt.close()

    ## Taylor diagrams (gradients)
    resfile=workpath+"/Taylor_diagram_grads_"+day+".png"
    label=['GT',\
           r"$\nabla_{OI} (nadir)$",r"$\nabla_{OI} (nadir+swot)$",\
           r"$\nabla_{AnDA}$ (nadir)",r"$\nabla_{AnDA} (nadir+swot)$",\
           r"$\nabla_{Post-AnDA}$ (nadir)",r"$\nabla_{Post-AnDA} (nadir+swot)$",\
           r"$\nabla_{FP-ConvAE}$ (nadir)",r"$\nabla_{FP-ConvAE} (nadir+swot)$",\
           r"$\nabla_{FP-GENN}$ (nadir)",r"$\nabla_{FP-GENN} (nadir+swot)$",\
           r"$\nabla_{VE-DINEOF}$ (nadir)",r"$\nabla_{VE-DINEOF} (nadir+swot)$"]
    series={'gt':gt,
            'Grad_OI_nadir': Grad_OI_nadir,'Grad_OI_nadirswot': Grad_OI_nadirswot,
            'Grad_AnDA_nadir':AnDA_nadir,'Grad_AnDA_nadirswot':AnDA_nadirswot,
            'Grad_Post_AnDA_nadir':Post_AnDA_nadir,'Grad_Post_AnDA_nadirswot':Post_AnDA_nadirswot,
            'Grad_FP_ConvAE_nadir':Grad_FP_ConvAE_nadir,'Grad_FP_ConvAE_nadirswot':Grad_FP_ConvAE_nadirswot,
            'Grad_FP_GENN_nadir':Grad_FP_GENN_nadir,'Grad_FP_GENN_nadirswot':Grad_FP_GENN_nadirswot,
            'Grad_VE_DINEOF_nadir':VE_DINEOF_nadir,'Grad_VE_DINEOF_nadirswot':VE_DINEOF_nadirswot}
    Taylor_diag(series,label,\
                styles=['k','s','p','o','p','o','p','o','p','o','p','o'],\
                colors=['k','y','mediumseagreen','mediumseagreen',\
                                'seagreen','seagreen',\
                                'steelblue','steelblue',\
                                'mediumorchid','mediumorchid',\
                                'darkorange','darkorange'])
    plt.savefig(resfile)
    plt.close()

    ## Radial Power Spectrum (RAPS)
    resfile=workpath+"/results_AnDA_RAPS_"+day+".png"
    f_ref, Pf_GT    			= raPsd2dv1(gt,resssh,True)
    f0_nadir, Pf_OI_nadir 		= raPsd2dv1(OI_nadir,resssh,True)
    f0_nadirswot, Pf_OI_nadirswot       = raPsd2dv1(OI_nadirswot,resssh,True)
    f1_nadir, Pf_AnDA_nadir  		= raPsd2dv1(AnDA_nadir,resssh,True)
    f1_nadirswot, Pf_AnDA_nadirswot  	= raPsd2dv1(AnDA_nadirswot,resssh,True)
    f2_nadir, Pf_postAnDA_nadir 	= raPsd2dv1(Post_AnDA_nadir,resssh,True)
    f2_nadirswot, Pf_postAnDA_nadirswot = raPsd2dv1(Post_AnDA_nadirswot,resssh,True)
    f3_nadir, Pf_VE_DINEOF_nadir   	= raPsd2dv1(VE_DINEOF_nadir,resssh,True)
    f3_nadirswot,Pf_VE_DINEOF_nadirswot = raPsd2dv1(VE_DINEOF_nadirswot,resssh,True)
    f4_nadir, Pf_FP_ConvAE_nadir         = raPsd2dv1(FP_ConvAE_nadir,resssh,True)
    f4_nadirswot, Pf_FP_ConvAE_nadirswot = raPsd2dv1(FP_ConvAE_nadirswot,resssh,True)
    f5_nadir, Pf_FP_GENN_nadir         = raPsd2dv1(FP_GENN_nadir,resssh,True)
    f5_nadirswot, Pf_FP_GENN_nadirswot = raPsd2dv1(FP_GENN_nadirswot,resssh,True)
    wf_ref	        = 1/f_ref
    wf0_nadir	        = 1/f0_nadir
    wf0_nadirswot       = 1/f0_nadirswot
    wf1_nadir        	= 1/f1_nadir
    wf1_nadirswot       = 1/f1_nadirswot
    wf2_nadir         	= 1/f2_nadir
    wf2_nadirswot       = 1/f2_nadirswot
    wf3_nadir         	= 1/f3_nadir
    wf3_nadirswot       = 1/f3_nadirswot
    wf4_nadir           = 1/f4_nadir
    wf4_nadirswot       = 1/f4_nadirswot
    wf5_nadir           = 1/f5_nadir
    wf5_nadirswot       = 1/f5_nadirswot
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(wf_ref[1:],Pf_GT[1:],label='GT')
    ax.plot(wf0_nadir[1:],Pf_OI_nadir[1:],linestyle='solid',color='red',linewidth=.5,label='OI (nadir)')
    ax.plot(wf0_nadirswot[1:],Pf_OI_nadirswot[1:],linestyle='dashdot',color='red',linewidth=.5,label='OI (nadir+swot)')
    ax.plot(wf1_nadir[1:],Pf_AnDA_nadir[1:],linestyle='solid',color='mediumseagreen',linewidth=.5,label='AnDA (nadir)')
    ax.plot(wf1_nadirswot[1:],Pf_AnDA_nadirswot[1:],linestyle='dashdot',color='mediumseagreen',linewidth=.5,label='AnDA (nadir+swot)')
    ax.plot(wf2_nadir[1:],Pf_postAnDA_nadir[1:],linestyle='solid',color='seagreen',linewidth=.5,label='post-AnDA (nadir)')
    ax.plot(wf2_nadirswot[1:],Pf_postAnDA_nadirswot[1:],linestyle='dashdot',color='seagreen',linewidth=.5,label='post-AnDA (nadir+swot)')
    ax.plot(wf3_nadir[1:],Pf_VE_DINEOF_nadir[1:],linestyle='solid',color='steelblue',linewidth=.5,label='VE-DINEOF (nadir)')
    ax.plot(wf3_nadirswot[1:],Pf_VE_DINEOF_nadirswot[1:],linestyle='dashdot',color='steelblue',linewidth=.5,label='VE-DINEOF (nadir+swot)')
    ax.plot(wf4_nadir[1:],Pf_FP_ConvAE_nadir[1:],linestyle='solid',color='mediumorchid',linewidth=.5,label='FP-ConvAE (nadir)')
    ax.plot(wf4_nadirswot[1:],Pf_FP_ConvAE_nadirswot[1:],linestyle='dashdot',color='mediumorchid',linewidth=.5,label='FP-ConvAE (nadir+swot)')
    ax.plot(wf5_nadir[1:],Pf_FP_GENN_nadir[1:],linestyle='solid',color='darkorange',linewidth=.5,label='FP-GENN (nadir)')
    ax.plot(wf5_nadirswot[1:],Pf_FP_GENN_nadirswot[1:],linestyle='dashdot',color='darkorange',linewidth=.5,label='FP-GENN (nadir+swot)')
    ax.set_xlabel("Wavenumber", fontweight='bold')
    ax.set_ylabel("Power spectral density (m2/(cy/km))", fontweight='bold')
    ax.set_xscale('log') ; ax.set_yscale('log')
    plt.legend(loc='best',prop=dict(size='small'),frameon=False)
    plt.xticks([50, 100, 200, 500, 1000], ["50km", "100km", "200km", "500km", "1000km"])
    ax.invert_xaxis()
    plt.grid(which='both', linestyle='--')
    plt.savefig(resfile)        # save the figure
    plt.close()         	# close the figure

    ## Radial Power Spectrum error (err_RAPS)
    resfile=workpath+"/results_diffAnDA_RAPS_"+day+".png"
    f0_nadir, Pf_OI_nadir               = err_raPsd2dv1(OI_nadir,gt,resssh,True)
    f0_nadirswot, Pf_OI_nadirswot       = err_raPsd2dv1(OI_nadirswot,gt,resssh,True)
    f1_nadir, Pf_AnDA_nadir             = err_raPsd2dv1(AnDA_nadir,gt,resssh,True)
    f1_nadirswot, Pf_AnDA_nadirswot     = err_raPsd2dv1(AnDA_nadirswot,gt,resssh,True)
    f2_nadir, Pf_postAnDA_nadir         = err_raPsd2dv1(Post_AnDA_nadir,gt,resssh,True)
    f2_nadirswot, Pf_postAnDA_nadirswot = err_raPsd2dv1(Post_AnDA_nadirswot,gt,resssh,True)
    f3_nadir, Pf_VE_DINEOF_nadir        = err_raPsd2dv1(VE_DINEOF_nadir,gt,resssh,True)
    f3_nadirswot,Pf_VE_DINEOF_nadirswot = err_raPsd2dv1(VE_DINEOF_nadirswot,gt,resssh,True)
    f4_nadir, Pf_FP_ConvAE_nadir         = err_raPsd2dv1(FP_ConvAE_nadir,gt,resssh,True)
    f4_nadirswot, Pf_FP_ConvAE_nadirswot = err_raPsd2dv1(FP_ConvAE_nadirswot,gt,resssh,True)
    f5_nadir, Pf_FP_GENN_nadir         = err_raPsd2dv1(FP_GENN_nadir,gt,resssh,True)
    f5_nadirswot, Pf_FP_GENN_nadirswot = err_raPsd2dv1(FP_GENN_nadirswot,gt,resssh,True)
    wf0_nadir           = 1/f0_nadir
    wf0_nadirswot       = 1/f0_nadirswot
    wf1_nadir           = 1/f1_nadir
    wf1_nadirswot       = 1/f1_nadirswot
    wf2_nadir           = 1/f2_nadir
    wf2_nadirswot       = 1/f2_nadirswot
    wf3_nadir           = 1/f3_nadir
    wf3_nadirswot       = 1/f3_nadirswot
    wf4_nadir           = 1/f4_nadir
    wf4_nadirswot       = 1/f4_nadirswot
    wf5_nadir           = 1/f5_nadir
    wf5_nadirswot       = 1/f5_nadirswot
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(wf0_nadir[1:],Pf_OI_nadir[1:],linestyle='solid',color='red',linewidth=.5,label='OI (nadir)')
    ax.plot(wf0_nadirswot[1:],Pf_OI_nadirswot[1:],linestyle='dashdot',color='red',linewidth=.5,label='OI (nadir+swot)')
    ax.plot(wf1_nadir[1:],Pf_AnDA_nadir[1:],linestyle='solid',color='mediumseagreen',linewidth=.5,label='AnDA (nadir)')
    ax.plot(wf1_nadirswot[1:],Pf_AnDA_nadirswot[1:],linestyle='dashdot',color='mediumseagreen',linewidth=.5,label='AnDA (nadir+swot)')
    ax.plot(wf2_nadir[1:],Pf_postAnDA_nadir[1:],linestyle='solid',color='seagreen',linewidth=.5,label='post-AnDA (nadir)')
    ax.plot(wf2_nadirswot[1:],Pf_postAnDA_nadirswot[1:],linestyle='dashdot',color='seagreen',linewidth=.5,label='post-AnDA (nadir+swot)')
    ax.plot(wf3_nadir[1:],Pf_VE_DINEOF_nadir[1:],linestyle='solid',color='steelblue',linewidth=.5,label='VE-DINEOF (nadir)')
    ax.plot(wf3_nadirswot[1:],Pf_VE_DINEOF_nadirswot[1:],linestyle='dashdot',color='steelblue',linewidth=.5,label='VE-DINEOF (nadir+swot)')
    ax.plot(wf4_nadir[1:],Pf_FP_ConvAE_nadir[1:],linestyle='solid',color='mediumorchid',linewidth=.5,label='FP-ConvAE (nadir)')
    ax.plot(wf4_nadirswot[1:],Pf_FP_ConvAE_nadirswot[1:],linestyle='dashdot',color='mediumorchid',linewidth=.5,label='FP-ConvAE (nadir+swot)')
    ax.plot(wf5_nadir[1:],Pf_FP_GENN_nadir[1:],linestyle='solid',color='darkorange',linewidth=.5,label='FP-GENN (nadir)')
    ax.plot(wf5_nadirswot[1:],Pf_FP_GENN_nadirswot[1:],linestyle='dashdot',color='darkorange',linewidth=.5,label='FP-GENN (nadir+swot)')
    ax.set_xlabel("Wavenumber", fontweight='bold')
    ax.set_ylabel("Power spectral density (m2/(cy/km))", fontweight='bold')
    ax.set_xscale('log') ; ax.set_yscale('log')
    plt.legend(loc='best',prop=dict(size='small'),frameon=False)
    plt.xticks([50, 100, 200, 500, 1000], ["50km", "100km", "200km", "500km", "1000km"])
    ax.invert_xaxis()
    plt.grid(which='both', linestyle='--')
    plt.savefig(resfile)        # save the figure
    plt.close()                 # close the figure

    ## Plot averaged RAPS (nadir)
    resfile=workpath+"/results_diffAnDA_avg_RAPS_nadir"+day+".png"
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(wf0_nadir[1:],Pf_OI_nadir[1:],linestyle='solid',color='red',linewidth=.5,label='OI (nadir)')
    ax.plot(wf1_nadir[1:],Pf_AnDA_nadir[1:],linestyle='solid',color='mediumseagreen',linewidth=.5,label='AnDA (nadir)')
    ax.plot(wf2_nadir[1:],Pf_postAnDA_nadir[1:],linestyle='solid',color='seagreen',linewidth=.5,label='post-AnDA (nadir)')
    ax.plot(wf3_nadir[1:],Pf_VE_DINEOF_nadir[1:],linestyle='solid',color='steelblue',linewidth=.5,label='VE-DINEOF (nadir)')
    ax.plot(wf4_nadir[1:],Pf_FP_ConvAE_nadir[1:],linestyle='solid',color='mediumorchid',linewidth=.5,label='FP-ConvAE (nadir)')
    ax.plot(wf5_nadir[1:],Pf_FP_GENN_nadir[1:],linestyle='solid',color='darkorange',linewidth=.5,label='FP-GENN (nadir)')
    ax.set_xlabel("Wavenumber", fontweight='bold')
    ax.set_ylabel("Power spectral density (m2/(cy/km))", fontweight='bold')
    ax.set_xscale('log') ; ax.set_yscale('log')
    plt.legend(loc='best',prop=dict(size='small'),frameon=False)
    plt.xticks([50, 100, 200, 500, 1000], ["50km", "100km", "200km", "500km", "1000km"])
    ax.invert_xaxis()
    plt.grid(which='both', linestyle='--')
    plt.savefig(resfile)# save the figure
    plt.close() # close the figure

    ## Plot averaged RAPS (nadirswot)
    resfile=workpath+"/results_diffAnDA_avg_RAPS_nadirswot"+day+".png"
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(wf0_nadirswot[1:],Pf_OI_nadirswot[1:],linestyle='solid',color='red',linewidth=.5,label='OI (nadir+swot)')
    ax.plot(wf1_nadirswot[1:],Pf_AnDA_nadirswot[1:],linestyle='solid',color='mediumseagreen',linewidth=.5,label='AnDA (nadirswot)')
    ax.plot(wf2_nadirswot[1:],Pf_postAnDA_nadirswot[1:],linestyle='solid',color='seagreen',linewidth=.5,label='post-AnDA (nadirswot)')
    ax.plot(wf3_nadirswot[1:],Pf_VE_DINEOF_nadirswot[1:],linestyle='solid',color='steelblue',linewidth=.5,label='VE-DINEOF (nadirswot)')
    ax.plot(wf4_nadirswot[1:],Pf_FP_ConvAE_nadirswot[1:],linestyle='solid',color='mediumorchid',linewidth=.5,label='FP-ConvAE (nadirswot)')
    ax.plot(wf5_nadirswot[1:],Pf_FP_GENN_nadirswot[1:],linestyle='solid',color='darkorange',linewidth=.5,label='FP-GENN (nadirswot)')
    ax.set_xlabel("Wavenumber", fontweight='bold')
    ax.set_ylabel("Power spectral density (m2/(cy/km))", fontweight='bold')
    ax.set_xscale('log') ; ax.set_yscale('log')
    plt.legend(loc='best',prop=dict(size='small'),frameon=False)
    plt.xticks([50, 100, 200, 500, 1000], ["50km", "100km", "200km", "500km", "1000km"])
    ax.invert_xaxis()
    plt.grid(which='both', linestyle='--')
    plt.savefig(resfile)# save the figure
    plt.close() # close the figure'''

## Taylor diagrams
resfile=workpath+"/Taylor_diagram.png"
label=['GT',\
       'OI (nadir)','OI (nadir+swot)',\
       'AnDA (nadir)','AnDA (nadir+swot)',\
       'Post-AnDA (nadir)','Post-AnDA (nadir+swot)',\
       'FP-ConvAE (nadir)','FP-ConvAE (nadir+swot)',\
       'FP-GENN (nadir)','FP-GENN (nadir+swot)',\
       'VE-DINEOF (nadir)','VE-DINEOF (nadir+swot)']
# apply HPF to visualize Taylor diagrams only for small scales
HR = AnDA_ssh_1.itrp_OI[:,:indLon,:indLat]
lr = np.copy(HR).reshape(HR.shape[0],-1)
tmp = lr[0,:]
sea_v2 = np.where(~np.isnan(tmp))[0]
lr_no_land = lr[:,sea_v2]
pca = PCA(n_components=.75)
score_global = pca.fit_transform(lr_no_land)
coeff_global = pca.components_.T
mu_global = pca.mean_
DataReconstructed_global = np.dot(score_global, coeff_global.T) +mu_global
lr[:,sea_v2] = DataReconstructed_global
lr = lr.reshape(HR.shape).flatten()
series={'gt':AnDA_ssh_1.GT[:,:indLon,:indLat].flatten()-lr,
        'OI_nadir':AnDA_ssh_1_nadir.itrp_OI[:,:indLon,:indLat].flatten()-lr,
        'OI_nadirswot':AnDA_ssh_1_nadirswot.itrp_OI[:,:indLon,:indLat].flatten()-lr,
        'AnDA_nadir':AnDA_ssh_1_nadir.itrp_AnDA[:,:indLon,:indLat].flatten()-lr,\
        'AnDA_nadirswot':AnDA_ssh_1_nadirswot.itrp_AnDA[:,:indLon,:indLat].flatten()-lr,
        'Post_AnDA_nadir':AnDA_ssh_1_nadir.itrp_postAnDA[:,:indLon,:indLat].flatten()-lr,\
        'Post_AnDA_nadirswot':AnDA_ssh_1_nadirswot.itrp_postAnDA[:,:indLon,:indLat].flatten()-lr,
        'FP_ConvAE_nadir':itrp_FP_ConvAE_nadir[:,:indLon,:indLat].flatten()-lr,\
        'FP_ConvAE_nadirswot':itrp_FP_ConvAE_nadirswot[:,:indLon,:indLat].flatten()-lr,
        'FP_GENN_nadir':itrp_FP_GENN_nadir[:,:indLon,:indLat].flatten()-lr,\
        'FP_GENN_nadirswot':itrp_FP_GENN_nadirswot[:,:indLon,:indLat].flatten()-lr,
        'VE_DINEOF_nadir':itrp_dineof_nadir[:,:indLon,:indLat].flatten()-lr,\
        'VE_DINEOF_nadirswot':itrp_dineof_nadirswot[:,:indLon,:indLat].flatten()-lr}
Taylor_diag(series,label,\
            styles=['k','s','p','o','p','o','p','o','p','o','p','o'],\
            colors=['k','y','mediumseagreen','mediumseagreen',\
                                'seagreen','seagreen',\
                                'steelblue','steelblue',\
                                'mediumorchid','mediumorchid',\
                                'darkorange','darkorange'])
plt.savefig(resfile)
plt.close()

## Plot averaged RAPS
resfile=workpath+"/results_AnDA_avg_RAPS.png"
f_ref, Pf_GT                         = avg_raPsd2dv1(AnDA_ssh_1.GT,resssh,True)
f0_nadir, Pf_OI_nadir                = avg_raPsd2dv1(AnDA_ssh_1_nadir.itrp_OI,resssh,True)
f0_nadirswot, Pf_OI_nadirswot        = avg_raPsd2dv1(AnDA_ssh_1_nadirswot.itrp_OI,resssh,True)
f1_nadir, Pf_AnDA_nadir              = avg_raPsd2dv1(AnDA_ssh_1_nadir.itrp_AnDA,resssh,True)
f1_nadirswot, Pf_AnDA_nadirswot      = avg_raPsd2dv1(AnDA_ssh_1_nadirswot.itrp_AnDA,resssh,True)
f2_nadir, Pf_postAnDA_nadir          = avg_raPsd2dv1(AnDA_ssh_1_nadir.itrp_postAnDA,resssh,True)
f2_nadirswot, Pf_postAnDA_nadirswot  = avg_raPsd2dv1(AnDA_ssh_1_nadirswot.itrp_postAnDA,resssh,True)
f3_nadir, Pf_VE_DINEOF_nadir         = avg_raPsd2dv1(itrp_dineof_nadir,resssh,True)
f3_nadirswot,Pf_VE_DINEOF_nadirswot  = avg_raPsd2dv1(itrp_dineof_nadirswot,resssh,True)
f4_nadir, Pf_FP_ConvAE_nadir         = avg_raPsd2dv1(itrp_FP_ConvAE_nadir,resssh,True)
f4_nadirswot, Pf_FP_ConvAE_nadirswot = avg_raPsd2dv1(itrp_FP_ConvAE_nadirswot,resssh,True)
f5_nadir, Pf_FP_GENN_nadir           = avg_raPsd2dv1(itrp_FP_GENN_nadir,resssh,True)
f5_nadirswot, Pf_FP_GENN_nadirswot   = avg_raPsd2dv1(itrp_FP_GENN_nadirswot,resssh,True)
wf_ref        = 1/f_ref
wf0_nadir     = 1/f0_nadir
wf0_nadirswot = 1/f0_nadirswot
wf1_nadir     = 1/f1_nadir
wf1_nadirswot = 1/f1_nadirswot
wf2_nadir     = 1/f2_nadir
wf2_nadirswot = 1/f2_nadirswot
wf3_nadir     = 1/f3_nadir
wf3_nadirswot = 1/f3_nadirswot
wf4_nadir     = 1/f4_nadir
wf4_nadirswot = 1/f4_nadirswot
wf5_nadir     = 1/f5_nadir
wf5_nadirswot = 1/f5_nadirswot
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(wf_ref[1:],Pf_GT[1:],label='GT')
ax.plot(wf0_nadir[1:],Pf_OI_nadir[1:],linestyle='solid',color='red',linewidth=.5,label='OI (nadir)')
ax.plot(wf0_nadirswot[1:],Pf_OI_nadirswot[1:],linestyle='dashdot',color='red',linewidth=.5,label='OI (nadir+swot)')
ax.plot(wf1_nadir[1:],Pf_AnDA_nadir[1:],linestyle='solid',color='mediumseagreen',linewidth=.5,label='AnDA (nadir)')
ax.plot(wf1_nadirswot[1:],Pf_AnDA_nadirswot[1:],linestyle='dashdot',color='mediumseagreen',linewidth=.5,label='AnDA (nadir+swot)')
ax.plot(wf2_nadir[1:],Pf_postAnDA_nadir[1:],linestyle='solid',color='seagreen',linewidth=.5,label='post-AnDA (nadir)')
ax.plot(wf2_nadirswot[1:],Pf_postAnDA_nadirswot[1:],linestyle='dashdot',color='seagreen',linewidth=.5,label='post-AnDA (nadir+swot)')
ax.plot(wf3_nadir[1:],Pf_VE_DINEOF_nadir[1:],linestyle='solid',color='steelblue',linewidth=.5,label='VE-DINEOF (nadir)')
ax.plot(wf3_nadirswot[1:],Pf_VE_DINEOF_nadirswot[1:],linestyle='dashdot',color='steelblue',linewidth=.5,label='VE-DINEOF (nadir+swot)')
ax.plot(wf4_nadir[1:],Pf_FP_ConvAE_nadir[1:],linestyle='solid',color='mediumorchid',linewidth=.5,label='FP-ConvAE (nadir)')
ax.plot(wf4_nadirswot[1:],Pf_FP_ConvAE_nadirswot[1:],linestyle='dashdot',color='mediumorchid',linewidth=.5,label='FP-ConvAE (nadir+swot)')
ax.plot(wf5_nadir[1:],Pf_FP_GENN_nadir[1:],linestyle='solid',color='darkorange',linewidth=.5,label='FP-GENN (nadir)')
ax.plot(wf5_nadirswot[1:],Pf_FP_GENN_nadirswot[1:],linestyle='dashdot',color='darkorange',linewidth=.5,label='FP-GENN (nadir+swot)')
ax.set_xlabel("Wavenumber", fontweight='bold')
ax.set_ylabel("Power spectral density (m2/(cy/km))", fontweight='bold')
ax.set_xscale('log') ; ax.set_yscale('log')
plt.legend(loc='best',prop=dict(size='small'),frameon=False)
plt.xticks([50, 100, 200, 500, 1000], ["50km", "100km", "200km", "500km", "1000km"])
ax.invert_xaxis()
plt.grid(which='both', linestyle='--')
plt.savefig(resfile)# save the figure
plt.close() # close the figure

## Plot averaged RAPS (nadir)
resfile=workpath+"/results_AnDA_avg_RAPS_nadir.png"
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(wf_ref[1:],Pf_GT[1:],label='GT')
ax.plot(wf0_nadir[1:],Pf_OI_nadir[1:],linestyle='solid',color='red',linewidth=.5,label='OI (nadir)')
ax.plot(wf1_nadir[1:],Pf_AnDA_nadir[1:],linestyle='solid',color='mediumseagreen',linewidth=.5,label='AnDA (nadir)')
ax.plot(wf2_nadir[1:],Pf_postAnDA_nadir[1:],linestyle='solid',color='seagreen',linewidth=.5,label='post-AnDA (nadir)')
ax.plot(wf3_nadir[1:],Pf_VE_DINEOF_nadir[1:],linestyle='solid',color='steelblue',linewidth=.5,label='VE-DINEOF (nadir)')
ax.plot(wf4_nadir[1:],Pf_FP_ConvAE_nadir[1:],linestyle='solid',color='mediumorchid',linewidth=.5,label='FP-ConvAE (nadir)')
ax.plot(wf5_nadir[1:],Pf_FP_GENN_nadir[1:],linestyle='solid',color='darkorange',linewidth=.5,label='FP-GENN (nadir)')
ax.set_xlabel("Wavenumber", fontweight='bold')
ax.set_ylabel("Power spectral density (m2/(cy/km))", fontweight='bold')
ax.set_xscale('log') ; ax.set_yscale('log')
plt.legend(loc='best',prop=dict(size='small'),frameon=False)
plt.xticks([50, 100, 200, 500, 1000], ["50km", "100km", "200km", "500km", "1000km"])
ax.invert_xaxis()
plt.grid(which='both', linestyle='--')
plt.savefig(resfile)# save the figure
plt.close() # close the figure

## Plot averaged RAPS (nadirswot)
resfile=workpath+"/results_AnDA_avg_RAPS_nadirswot.png"
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(wf_ref[1:],Pf_GT[1:],label='GT')
ax.plot(wf0_nadirswot[1:],Pf_OI_nadirswot[1:],linestyle='solid',color='red',linewidth=.5,label='OI (nadir+swot)')
ax.plot(wf1_nadirswot[1:],Pf_AnDA_nadirswot[1:],linestyle='solid',color='mediumseagreen',linewidth=.5,label='AnDA (nadir+swot)')
ax.plot(wf2_nadirswot[1:],Pf_postAnDA_nadirswot[1:],linestyle='solid',color='seagreen',linewidth=.5,label='post-AnDA (nadir+swot)')
ax.plot(wf3_nadirswot[1:],Pf_VE_DINEOF_nadirswot[1:],linestyle='solid',color='steelblue',linewidth=.5,label='VE-DINEOF (nadir+swot)')
ax.plot(wf4_nadirswot[1:],Pf_FP_ConvAE_nadirswot[1:],linestyle='solid',color='mediumorchid',linewidth=.5,label='FP-ConvAE (nadir+swot)')
ax.plot(wf5_nadirswot[1:],Pf_FP_GENN_nadirswot[1:],linestyle='solid',color='darkorange',linewidth=.5,label='FP-GENN (nadir+swot)')
ax.set_xlabel("Wavenumber", fontweight='bold')
ax.set_ylabel("Power spectral density (m2/(cy/km))", fontweight='bold')
ax.set_xscale('log') ; ax.set_yscale('log')
plt.legend(loc='best',prop=dict(size='small'),frameon=False)
plt.xticks([50, 100, 200, 500, 1000], ["50km", "100km", "200km", "500km", "1000km"])
ax.invert_xaxis()
plt.grid(which='both', linestyle='--')
plt.savefig(resfile)# save the figure
plt.close() # close the figure

## Plot averaged normalize error RAPS
resfile=workpath+"/results_diffAnDA_avg_RAPS.png"
f0_nadir, Pf_OI_nadir                = avg_err_raPsd2dv1(AnDA_ssh_1_nadir.itrp_OI[:,:indLon,:indLat],AnDA_ssh_1.GT[:,:indLon,:indLat],resssh,True)
f0_nadirswot, Pf_OI_nadirswot        = avg_err_raPsd2dv1(AnDA_ssh_1_nadirswot.itrp_OI[:,:indLon,:indLat],AnDA_ssh_1.GT[:,:indLon,:indLat],gt,resssh,True)
f1_nadir, Pf_AnDA_nadir              = avg_err_raPsd2dv1(AnDA_ssh_1_nadir.itrp_AnDA[:,:indLon,:indLat],AnDA_ssh_1.GT[:,:indLon,:indLat],resssh,True)
f1_nadirswot, Pf_AnDA_nadirswot      = avg_err_raPsd2dv1(AnDA_ssh_1_nadirswot.itrp_AnDA[:,:indLon,:indLat],AnDA_ssh_1.GT[:,:indLon,:indLat],resssh,True)
f2_nadir, Pf_postAnDA_nadir          = avg_err_raPsd2dv1(AnDA_ssh_1_nadir.itrp_postAnDA[:,:indLon,:indLat],AnDA_ssh_1.GT[:,:indLon,:indLat],resssh,True)
f2_nadirswot, Pf_postAnDA_nadirswot  = avg_err_raPsd2dv1(AnDA_ssh_1_nadirswot.itrp_postAnDA[:,:indLon,:indLat],AnDA_ssh_1.GT[:,:indLon,:indLat],resssh,True)
f3_nadir, Pf_VE_DINEOF_nadir         = avg_err_raPsd2dv1(itrp_dineof_nadir[:,:indLon,:indLat],AnDA_ssh_1.GT[:,:indLon,:indLat],resssh,True)
f3_nadirswot,Pf_VE_DINEOF_nadirswot  = avg_err_raPsd2dv1(itrp_dineof_nadirswot[:,:indLon,:indLat],AnDA_ssh_1.GT[:,:indLon,:indLat],resssh,True)
f4_nadir, Pf_FP_ConvAE_nadir         = avg_err_raPsd2dv1(itrp_FP_ConvAE_nadir[:,:indLon,:indLat],AnDA_ssh_1.GT[:,:indLon,:indLat],resssh,True)
f4_nadirswot, Pf_FP_ConvAE_nadirswot = avg_err_raPsd2dv1(itrp_FP_ConvAE_nadirswot[:,:indLon,:indLat],AnDA_ssh_1.GT[:,:indLon,:indLat],resssh,True)
f5_nadir, Pf_FP_GENN_nadir           = avg_err_raPsd2dv1(itrp_FP_GENN_nadir[:,:indLon,:indLat],AnDA_ssh_1.GT[:,:indLon,:indLat],resssh,True)
f5_nadirswot, Pf_FP_GENN_nadirswot   = avg_err_raPsd2dv1(itrp_FP_GENN_nadirswot[:,:indLon,:indLat],AnDA_ssh_1.GT[:,:indLon,:indLat],resssh,True)
wf0_nadir = 1/f0_nadir
wf0_nadirswot = 1/f0_nadirswot
wf1_nadir   = 1/f1_nadir
wf1_nadirswot   = 1/f1_nadirswot
wf2_nadir   = 1/f2_nadir
wf2_nadirswot   = 1/f2_nadirswot
wf3_nadir   = 1/f3_nadir
wf3_nadirswot   = 1/f3_nadirswot
wf4_nadir           = 1/f4_nadir
wf4_nadirswot       = 1/f4_nadirswot
wf5_nadir           = 1/f5_nadir
wf5_nadirswot       = 1/f5_nadirswot
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(wf0_nadir[1:],Pf_OI_nadir[1:],linestyle='solid',color='red',linewidth=.5,label='OI (nadir)')
ax.plot(wf0_nadirswot[1:],Pf_OI_nadirswot[1:],linestyle='dashdot',color='red',linewidth=.5,label='OI (nadir+swot)')
ax.plot(wf1_nadir[1:],Pf_AnDA_nadir[1:],linestyle='solid',color='mediumseagreen',linewidth=.5,label='AnDA (nadir)')
ax.plot(wf1_nadirswot[1:],Pf_AnDA_nadirswot[1:],linestyle='dashdot',color='mediumseagreen',linewidth=.5,label='AnDA (nadir+swot)')
ax.plot(wf2_nadir[1:],Pf_postAnDA_nadir[1:],linestyle='solid',color='seagreen',linewidth=.5,label='post-AnDA (nadir)')
ax.plot(wf2_nadirswot[1:],Pf_postAnDA_nadirswot[1:],linestyle='dashdot',color='seagreen',linewidth=.5,label='post-AnDA (nadir+swot)')
ax.plot(wf3_nadir[1:],Pf_VE_DINEOF_nadir[1:],linestyle='solid',color='steelblue',linewidth=.5,label='VE-DINEOF (nadir)')
ax.plot(wf3_nadirswot[1:],Pf_VE_DINEOF_nadirswot[1:],linestyle='dashdot',color='steelblue',linewidth=.5,label='VE-DINEOF (nadir+swot)')
ax.plot(wf4_nadir[1:],Pf_FP_ConvAE_nadir[1:],linestyle='solid',color='mediumorchid',linewidth=.5,label='FP-ConvAE (nadir)')
ax.plot(wf4_nadirswot[1:],Pf_FP_ConvAE_nadirswot[1:],linestyle='dashdot',color='mediumorchid',linewidth=.5,label='FP-ConvAE (nadir+swot)')
ax.plot(wf5_nadir[1:],Pf_FP_GENN_nadir[1:],linestyle='solid',color='darkorange',linewidth=.5,label='FP-GENN (nadir)')
ax.plot(wf5_nadirswot[1:],Pf_FP_GENN_nadirswot[1:],linestyle='dashdot',color='darkorange',linewidth=.5,label='FP-GENN (nadir+swot)')
ax.set_xlabel("Wavenumber", fontweight='bold')
ax.set_ylabel("Power spectral density (m2/(cy/km))", fontweight='bold')
ax.set_xscale('log') ; ax.set_yscale('log')
plt.legend(loc='best',prop=dict(size='small'),frameon=False)
plt.xticks([50, 100, 200, 500, 1000], ["50km", "100km", "200km", "500km", "1000km"])
ax.invert_xaxis()
plt.grid(which='both', linestyle='--')
plt.savefig(resfile)# save the figure
plt.close() # close the figure

## Plot averaged normalize error RAPS (nadir)
resfile=workpath+"/results_diffAnDA_avg_RAPS_nadir.png"
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(wf0_nadir[1:],Pf_OI_nadir[1:],linestyle='solid',color='red',linewidth=.5,label='OI (nadir)')
ax.plot(wf1_nadir[1:],Pf_AnDA_nadir[1:],linestyle='solid',color='mediumseagreen',linewidth=.5,label='AnDA (nadir)')
ax.plot(wf2_nadir[1:],Pf_postAnDA_nadir[1:],linestyle='solid',color='seagreen',linewidth=.5,label='post-AnDA (nadir)')
ax.plot(wf3_nadir[1:],Pf_VE_DINEOF_nadir[1:],linestyle='solid',color='steelblue',linewidth=.5,label='VE-DINEOF (nadir)')
ax.plot(wf4_nadir[1:],Pf_FP_ConvAE_nadir[1:],linestyle='solid',color='mediumorchid',linewidth=.5,label='FP-ConvAE (nadir)')
ax.plot(wf5_nadir[1:],Pf_FP_GENN_nadir[1:],linestyle='solid',color='darkorange',linewidth=.5,label='FP-GENN (nadir)')
ax.set_xlabel("Wavenumber", fontweight='bold')
ax.set_ylabel("Power spectral density (m2/(cy/km))", fontweight='bold')
ax.set_xscale('log') ; ax.set_yscale('log')
plt.legend(loc='best',prop=dict(size='small'),frameon=False)
plt.xticks([50, 100, 200, 500, 1000], ["50km", "100km", "200km", "500km", "1000km"])
ax.invert_xaxis()
plt.grid(which='both', linestyle='--')
plt.savefig(resfile)# save the figure
plt.close() # close the figure


## Plot averaged normalize error RAPS (nadirswot)
resfile=workpath+"/results_diffAnDA_avg_RAPS_nadirswot.png"
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(wf0_nadirswot[1:],Pf_OI_nadirswot[1:],linestyle='solid',color='red',linewidth=.5,label='OI (nadir+swot)')
ax.plot(wf1_nadirswot[1:],Pf_AnDA_nadirswot[1:],linestyle='solid',color='mediumseagreen',linewidth=.5,label='AnDA (nadir+swot)')
ax.plot(wf2_nadirswot[1:],Pf_postAnDA_nadirswot[1:],linestyle='solid',color='seagreen',linewidth=.5,label='post-AnDA (nadir+swot)')
ax.plot(wf3_nadirswot[1:],Pf_VE_DINEOF_nadirswot[1:],linestyle='solid',color='steelblue',linewidth=.5,label='VE-DINEOF (nadir+swot)')
ax.plot(wf4_nadirswot[1:],Pf_FP_ConvAE_nadirswot[1:],linestyle='solid',color='mediumorchid',linewidth=.5,label='FP-ConvAE (nadir+swot)')
ax.plot(wf5_nadirswot[1:],Pf_FP_GENN_nadirswot[1:],linestyle='solid',color='darkorange',linewidth=.5,label='FP-GENN (nadir+swot)')
ax.set_xlabel("Wavenumber", fontweight='bold')
ax.set_ylabel("Power spectral density (m2/(cy/km))", fontweight='bold')
ax.set_xscale('log') ; ax.set_yscale('log')
plt.legend(loc='best',prop=dict(size='small'),frameon=False)
plt.xticks([50, 100, 200, 500, 1000], ["50km", "100km", "200km", "500km", "1000km"])
ax.invert_xaxis()
plt.grid(which='both', linestyle='--')
plt.savefig(resfile)# save the figure
plt.close() # close the figure

## Plot time series (all) 
ymax_ = np.ceil(np.max([nrmse_OI_nadir, nrmse_OI_nadirswot, nrmse_Post_AnDA_nadir, nrmse_VE_DINEOF_nadir, nrmse_FP_ConvAE_nadir, nrmse_FP_GENN_nadir,\
                nrmse_Post_AnDA_nadirswot,nrmse_VE_DINEOF_nadirswot,nrmse_FP_ConvAE_nadirswot,nrmse_FP_GENN_nadirswot])*100)/100

N = len(lday)
print(N)
# first axis with nRMSE time series
plt.plot(range(N),nrmse_OI_nadir,linestyle='solid',color='red',linewidth=2,label='OI (nadir)')
plt.plot(range(N),nrmse_Post_AnDA_nadir,linestyle='solid',color='seagreen',linewidth=1,label='post-AnDA (nadir)')
plt.plot(range(N),nrmse_VE_DINEOF_nadir,linestyle='solid',color='steelblue',linewidth=1,markerSize=2,label='VE-DINEOF (nadir)')
plt.plot(range(N),nrmse_FP_ConvAE_nadir,linestyle='solid',color='mediumorchid',linewidth=1,markerSize=2,label='FP-ConvAE (nadir)')
plt.plot(range(N),nrmse_FP_GENN_nadir,linestyle='solid',color='darkorange',linewidth=1,markerSize=2,label='FP-GENN (nadir)')
plt.plot(range(N),nrmse_OI_nadirswot,linestyle='dashdot',color='red',linewidth=2,label='OI (nadir+swot)')
plt.plot(range(N),nrmse_Post_AnDA_nadirswot,linestyle='dashdot',color='seagreen',linewidth=1,label='post-AnDA (nadir+swot)')
plt.plot(range(N),nrmse_VE_DINEOF_nadirswot,linestyle='dashdot',color='steelblue',linewidth=1,markerSize=2,label='VE-DINEOF (nadir+swot)')
plt.plot(range(N),nrmse_FP_ConvAE_nadirswot,linestyle='dashdot',color='mediumorchid',linewidth=1,markerSize=2,label='FP-ConvAE (nadir+swot)')
plt.plot(range(N),nrmse_FP_GENN_nadirswot,linestyle='dashdot',color='darkorange',linewidth=1,markerSize=2,label='FP-GENN (nadir+swot)')
# add vertical bar to divide the 4 periods
plt.axvline(x=19)
plt.axvline(x=39)
plt.axvline(x=59)
# graphical options
plt.ylim(0,ymax_)
plt.ylabel('nRMSE')
plt.xlabel('Time (days)')
plt.xticks([0,20,40,60],\
           [lday[0],lday[20],lday[40],lday[60]],\
           rotation=45, ha='right')
plt.margins(x=0)
plt.grid(True,alpha=.3)
plt.legend(loc='upper left',prop=dict(size='small'),frameon=False,bbox_to_anchor=(0,1.02,1,0.2),ncol=3,mode="expand")
# second axis with spatial coverage
axes2 = plt.twinx()
width=0.75
p1 = axes2.bar(range(N), nadcov, width,color='r',alpha=0.25)
p2 = axes2.bar(range(N), swotcov, width,bottom=nadcov,color='g',alpha=0.25)
axes2.set_ylim(0, 1)
axes2.set_ylabel('Spatial Coverage (%)')
axes2.margins(x=0)
resfile=workpath+"/TS_AnDA_nRMSE.png"
plt.savefig(resfile,bbox_inches="tight")    # save the figure
plt.close()         	# close the figure

## Plot time series (nadir)
N = len(lday)
print(N)
# first axis with nRMSE time series
plt.plot(range(N),nrmse_OI_nadir,linestyle='solid',color='red',linewidth=2,label='OI (nadir)')
plt.plot(range(N),nrmse_Post_AnDA_nadir,linestyle='solid',color='seagreen',linewidth=1,label='post-AnDA (nadir)')
plt.plot(range(N),nrmse_VE_DINEOF_nadir,linestyle='solid',color='steelblue',linewidth=1,markerSize=2,label='VE-DINEOF (nadir)')
plt.plot(range(N),nrmse_FP_ConvAE_nadir,linestyle='solid',color='mediumorchid',linewidth=1,markerSize=2,label='FP-ConvAE (nadir)')
plt.plot(range(N),nrmse_FP_GENN_nadir,linestyle='solid',color='darkorange',linewidth=1,markerSize=2,label='FP-GENN (nadir)')
# add vertical bar to divide the 4 periods
plt.axvline(x=19)
plt.axvline(x=39)
plt.axvline(x=59)
# graphical options
plt.ylim(0,ymax_)
plt.ylabel('nRMSE')
plt.xlabel('Time (days)')
plt.xticks([0,20,40,60],\
           [lday[0],lday[20],lday[40],lday[60]],\
           rotation=45, ha='right')
plt.margins(x=0)
plt.grid(True,alpha=.3)
plt.legend(loc='upper left',prop=dict(size='small'),frameon=False,bbox_to_anchor=(0,1.02,1,0.2),ncol=3,mode="expand")
# second axis with spatial coverage
axes2 = plt.twinx()
width=0.75
p1 = axes2.bar(range(N), nadcov, width,color='r',alpha=0.25)
p2 = axes2.bar(range(N), swotcov, width,bottom=nadcov,color='g',alpha=0.25)
axes2.set_ylim(0, 1)
axes2.set_ylabel('Spatial Coverage (%)')
axes2.margins(x=0)
resfile=workpath+"/TS_AnDA_nRMSE_nadir.png"
plt.savefig(resfile,bbox_inches="tight")    # save the figure
plt.close()             # close the figure

## Plot time series (nadir/swot)
N = len(lday)
print(N)
# first axis with nRMSE time series
plt.plot(range(N),nrmse_OI_nadirswot,linestyle='solid',color='red',linewidth=2,label='OI (nadir+swot)')
plt.plot(range(N),nrmse_Post_AnDA_nadirswot,linestyle='solid',color='seagreen',linewidth=1,label='post-AnDA (nadir+swot)')
plt.plot(range(N),nrmse_VE_DINEOF_nadirswot,linestyle='solid',color='steelblue',linewidth=1,markerSize=2,label='VE-DINEOF (nadir+swot)')
plt.plot(range(N),nrmse_FP_ConvAE_nadirswot,linestyle='solid',color='mediumorchid',linewidth=1,markerSize=2,label='FP-ConvAE (nadir+swot)')
plt.plot(range(N),nrmse_FP_GENN_nadirswot,linestyle='solid',color='darkorange',linewidth=1,markerSize=2,label='FP-GENN (nadir+swot)')
# add vertical bar to divide the 4 periods
plt.axvline(x=19)
plt.axvline(x=39)
plt.axvline(x=59)
# graphical options
plt.ylim(0,ymax_)
plt.ylabel('nRMSE')
plt.xlabel('Time (days)')
plt.xticks([0,20,40,60],\
           [lday[0],lday[20],lday[40],lday[60]],\
           rotation=45, ha='right')
plt.margins(x=0)
plt.grid(True,alpha=.3)
plt.legend(loc='upper left',prop=dict(size='small'),frameon=False,bbox_to_anchor=(0,1.02,1,0.2),ncol=3,mode="expand")
# second axis with spatial coverage
axes2 = plt.twinx()
width=0.75
p1 = axes2.bar(range(N), nadcov, width,color='r',alpha=0.25)
p2 = axes2.bar(range(N), swotcov, width,bottom=nadcov,color='g',alpha=0.25)
axes2.set_ylim(0, 1)
axes2.set_ylabel('Spatial Coverage (%)')
axes2.margins(x=0)
resfile=workpath+"/TS_AnDA_nRMSE_nadirswot.png"
plt.savefig(resfile,bbox_inches="tight")    # save the figure
plt.close()             # close the figure



