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
domain   = sys.argv[4] 
workpath = "/home3/scratch/mbeaucha/"+domain+"/scores_allmethods_AnDAnadlag_"+AnDA_lag+"_NNnadlag_"+NN_lag+"_"+type_obs
scratchpath = '/home3/scratch/mbeaucha/'+domain
if not os.path.exists(workpath):
    mk_dir_recursive(workpath)
#else:
#    shutil.rmtree(workpath)
#    mk_dir_recursive(workpath)    

# Reload saved AnDA result
file_results_nadir=scratchpath+'/resAnDA_nadir_nadlag_'+AnDA_lag+"_"+type_obs+'/saved_path.pickle'
with open(file_results_nadir, 'rb') as handle:
    AnDA_ssh_1, itrp_dineof = pickle.load(handle)
    AnDA_ssh_1_nadir = AnDA_ssh_1  
    itrp_dineof_nadir = itrp_dineof
file_results_nadirswot=scratchpath+'/resAnDA_nadirswot_nadlag_'+AnDA_lag+"_"+type_obs+'/saved_path.pickle'
with open(file_results_nadirswot, 'rb') as handle:
    AnDA_ssh_1, itrp_dineof = pickle.load(handle)
    AnDA_ssh_1_nadirswot = AnDA_ssh_1  
    itrp_dineof_nadirswot = itrp_dineof

# Reload saved ConvAE and GE-NN results
file_results_nadir=scratchpath+'/resIA_nadir_nadlag_'+NN_lag+"_"+type_obs+'/FP_ConvAE_womissing_wocov/saved_path_019_FP_ConvAE_womissing.pickle'
with open(file_results_nadir, 'rb') as handle:
    itrp_FP_ConvAE_nadir, rec_FP_ConvAE_nadir = pickle.load(handle)[2:]
file_results_nadirswot=scratchpath+'/resIA_nadirswot_nadlag_'+NN_lag+"_"+type_obs+'/FP_ConvAE_womissing_wocov/saved_path_019_FP_ConvAE_womissing.pickle'
with open(file_results_nadirswot, 'rb') as handle:
    itrp_FP_ConvAE_nadirswot, rec_FP_ConvAE_nadirswot = pickle.load(handle)[2:]
file_results_nadir=scratchpath+'/resIA_nadir_nadlag_'+NN_lag+"_"+type_obs+'/FP_GENN_wmissing_wOI/saved_path_019_FP_GENN_wmissing.pickle'
with open(file_results_nadir, 'rb') as handle:
    itrp_FP_GENN_nadir, rec_FP_GENN_nadir = pickle.load(handle)[2:]
file_results_nadirswot=scratchpath+'/resIA_nadirswot_nadlag_'+NN_lag+"_"+type_obs+'/FP_GENN_wwmissing_wocov/saved_path_019_FP_GENN_wwmissing.pickle'
with open(file_results_nadirswot, 'rb') as handle:
    itrp_FP_GENN_nadirswot, rec_FP_GENN_nadirswot = pickle.load(handle)[2:]

			#*****************#
			# Display results #
			#*****************#
resssh = 4
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

##*** INIT SSH STATISTICS ***##
## Init variables for temporal analysis (R scores)
nrmse_OI_nadir=np.zeros(len(AnDA_ssh_1.GT))
nrmse_OI_nadirswot=np.zeros(len(AnDA_ssh_1.GT))
nrmse_Post_AnDA_nadir=np.zeros(len(AnDA_ssh_1.GT))
nrmse_VE_DINEOF_nadir=np.zeros(len(AnDA_ssh_1.GT))
nrmse_FP_ConvAE_nadir=np.zeros(len(AnDA_ssh_1.GT))
nrmse_FP_GENN_nadir=np.zeros(len(AnDA_ssh_1.GT))
nrmse_Post_AnDA_nadirswot=np.zeros(len(AnDA_ssh_1.GT))
nrmse_VE_DINEOF_nadirswot=np.zeros(len(AnDA_ssh_1.GT))
nrmse_FP_ConvAE_nadirswot=np.zeros(len(AnDA_ssh_1.GT))
nrmse_FP_GENN_nadirswot=np.zeros(len(AnDA_ssh_1.GT))
## Init variables for temporal analysis (R scores)
R_OI_nadir=np.zeros(len(AnDA_ssh_1.GT))
R_OI_nadirswot=np.zeros(len(AnDA_ssh_1.GT))
R_Post_AnDA_nadir=np.zeros(len(AnDA_ssh_1.GT))
R_VE_DINEOF_nadir=np.zeros(len(AnDA_ssh_1.GT))
R_FP_ConvAE_nadir=np.zeros(len(AnDA_ssh_1.GT))
R_FP_GENN_nadir=np.zeros(len(AnDA_ssh_1.GT))
R_Post_AnDA_nadirswot=np.zeros(len(AnDA_ssh_1.GT))
R_VE_DINEOF_nadirswot=np.zeros(len(AnDA_ssh_1.GT))
R_FP_ConvAE_nadirswot=np.zeros(len(AnDA_ssh_1.GT))
R_FP_GENN_nadirswot=np.zeros(len(AnDA_ssh_1.GT))
## Init variables for temporal analysis (I scores)
I_OI_nadir=np.zeros(len(AnDA_ssh_1.GT))
I_OI_nadirswot=np.zeros(len(AnDA_ssh_1.GT))
I_Post_AnDA_nadir=np.zeros(len(AnDA_ssh_1.GT))
I_VE_DINEOF_nadir=np.zeros(len(AnDA_ssh_1.GT))
I_FP_ConvAE_nadir=np.zeros(len(AnDA_ssh_1.GT))
I_FP_GENN_nadir=np.zeros(len(AnDA_ssh_1.GT))
I_Post_AnDA_nadirswot=np.zeros(len(AnDA_ssh_1.GT))
I_VE_DINEOF_nadirswot=np.zeros(len(AnDA_ssh_1.GT))
I_FP_ConvAE_nadirswot=np.zeros(len(AnDA_ssh_1.GT))
I_FP_GENN_nadirswot=np.zeros(len(AnDA_ssh_1.GT))
## Init variables for temporal analysis (AE scores)
AE_FP_ConvAE_nadir=np.zeros(len(AnDA_ssh_1.GT))
AE_FP_GENN_nadir=np.zeros(len(AnDA_ssh_1.GT))
AE_FP_ConvAE_nadirswot=np.zeros(len(AnDA_ssh_1.GT))
AE_FP_GENN_nadirswot=np.zeros(len(AnDA_ssh_1.GT))

##*** INIT GradSSH STATISTICS ***##
## Init variables for temporal analysis (R scores)
nrmse_Grad_OI_nadir=np.zeros(len(AnDA_ssh_1.GT))
nrmse_Grad_OI_nadirswot=np.zeros(len(AnDA_ssh_1.GT))
nrmse_Grad_Post_AnDA_nadir=np.zeros(len(AnDA_ssh_1.GT))
nrmse_Grad_VE_DINEOF_nadir=np.zeros(len(AnDA_ssh_1.GT))
nrmse_Grad_FP_ConvAE_nadir=np.zeros(len(AnDA_ssh_1.GT))
nrmse_Grad_FP_GENN_nadir=np.zeros(len(AnDA_ssh_1.GT))
nrmse_Grad_Post_AnDA_nadirswot=np.zeros(len(AnDA_ssh_1.GT))
nrmse_Grad_VE_DINEOF_nadirswot=np.zeros(len(AnDA_ssh_1.GT))
nrmse_Grad_FP_ConvAE_nadirswot=np.zeros(len(AnDA_ssh_1.GT))
nrmse_Grad_FP_GENN_nadirswot=np.zeros(len(AnDA_ssh_1.GT))
## Init variables for temporal analysis (R scores)
R_Grad_OI_nadir=np.zeros(len(AnDA_ssh_1.GT))
R_Grad_OI_nadirswot=np.zeros(len(AnDA_ssh_1.GT))
R_Grad_Post_AnDA_nadir=np.zeros(len(AnDA_ssh_1.GT))
R_Grad_VE_DINEOF_nadir=np.zeros(len(AnDA_ssh_1.GT))
R_Grad_FP_ConvAE_nadir=np.zeros(len(AnDA_ssh_1.GT))
R_Grad_FP_GENN_nadir=np.zeros(len(AnDA_ssh_1.GT))
R_Grad_Post_AnDA_nadirswot=np.zeros(len(AnDA_ssh_1.GT))
R_Grad_VE_DINEOF_nadirswot=np.zeros(len(AnDA_ssh_1.GT))
R_Grad_FP_ConvAE_nadirswot=np.zeros(len(AnDA_ssh_1.GT))
R_Grad_FP_GENN_nadirswot=np.zeros(len(AnDA_ssh_1.GT))
## Init variables for temporal analysis (I scores)
I_Grad_OI_nadir=np.zeros(len(AnDA_ssh_1.GT))
I_Grad_OI_nadirswot=np.zeros(len(AnDA_ssh_1.GT))
I_Grad_Post_AnDA_nadir=np.zeros(len(AnDA_ssh_1.GT))
I_Grad_VE_DINEOF_nadir=np.zeros(len(AnDA_ssh_1.GT))
I_Grad_FP_ConvAE_nadir=np.zeros(len(AnDA_ssh_1.GT))
I_Grad_FP_GENN_nadir=np.zeros(len(AnDA_ssh_1.GT))
I_Grad_Post_AnDA_nadirswot=np.zeros(len(AnDA_ssh_1.GT))
I_Grad_VE_DINEOF_nadirswot=np.zeros(len(AnDA_ssh_1.GT))
I_Grad_FP_ConvAE_nadirswot=np.zeros(len(AnDA_ssh_1.GT))
I_Grad_FP_GENN_nadirswot=np.zeros(len(AnDA_ssh_1.GT))
## Init variables for temporal analysis (AE scores)
AE_Grad_FP_ConvAE_nadir=np.zeros(len(AnDA_ssh_1.GT))
AE_Grad_FP_GENN_nadir=np.zeros(len(AnDA_ssh_1.GT))
AE_Grad_FP_ConvAE_nadirswot=np.zeros(len(AnDA_ssh_1.GT))
AE_Grad_FP_GENN_nadirswot=np.zeros(len(AnDA_ssh_1.GT))

## Observations spatial coverage
nadcov=np.zeros(len(AnDA_ssh_1.GT))
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
    ## Maps
    resfile1=workpath+"/results_"+day+".png"
    resfile2=workpath+"/results_Grad_"+day+".png"
    # Load data
    gt 				= AnDA_ssh_1.GT[i,:indLat,:indLon]
    Grad_gt             	= Gradient(gt,2)
    # nadir
    OI_nadir                    = AnDA_ssh_1_nadir.itrp_OI[i,:indLat,:indLon]
    Grad_OI_nadir               = Gradient(OI_nadir,2)
    obs_nadir 			= AnDA_ssh_1_nadir.Obs[i,:indLat,:indLon]
    mask1_nadir                 = np.where(np.isnan(obs_nadir),np.nan,1)
    mask2_nadir                 = np.where(np.isnan(obs_nadir),1,np.nan)
    VE_DINEOF_nadir           	= itrp_dineof_nadir[i,:indLat,:indLon]
    Grad_VE_DINEOF_nadir 	= Gradient(VE_DINEOF_nadir,2)
    AnDA_nadir 		 	= AnDA_ssh_1_nadir.itrp_AnDA[i,:indLat,:indLon]
    Grad_AnDA_nadir          	= Gradient(AnDA_nadir,2)
    Post_AnDA_nadir 		= AnDA_ssh_1_nadir.itrp_postAnDA[i,:indLat,:indLon]
    Grad_Post_AnDA_nadir	= Gradient(Post_AnDA_nadir,2)
    FP_ConvAE_nadir		= itrp_FP_ConvAE_nadir[i,:indLat,:indLon]
    Grad_FP_ConvAE_nadir        = Gradient(FP_ConvAE_nadir,2)
    FP_GENN_nadir               = itrp_FP_GENN_nadir[i,:indLat,:indLon]
    Grad_FP_GENN_nadir          = Gradient(FP_GENN_nadir,2)
    rFP_ConvAE_nadir            = rec_FP_ConvAE_nadir[i,:indLat,:indLon]
    rGrad_FP_ConvAE_nadir       = Gradient(rFP_ConvAE_nadir,2)
    rFP_GENN_nadir              = rec_FP_GENN_nadir[i,:indLat,:indLon]
    rGrad_FP_GENN_nadir         = Gradient(rFP_GENN_nadir,2)
    # nadirswot
    OI_nadirswot                = AnDA_ssh_1_nadirswot.itrp_OI[i,:indLat,:indLon]
    Grad_OI_nadirswot           = Gradient(OI_nadirswot,2)
    obs_nadirswot 		= AnDA_ssh_1_nadirswot.Obs[i,:indLat,:indLon]
    mask1_nadirswot             = np.where(np.isnan(obs_nadirswot),np.nan,1)
    mask2_nadirswot             = np.where(np.isnan(obs_nadirswot),1,np.nan) 
    VE_DINEOF_nadirswot         = itrp_dineof_nadirswot[i,:indLat,:indLon]
    Grad_VE_DINEOF_nadirswot    = Gradient(VE_DINEOF_nadirswot,2)
    AnDA_nadirswot 		= AnDA_ssh_1_nadirswot.itrp_AnDA[i,:indLat,:indLon]
    Grad_AnDA_nadirswot         = Gradient(AnDA_nadirswot,2)
    Post_AnDA_nadirswot         = AnDA_ssh_1_nadirswot.itrp_postAnDA[i,:indLat,:indLon]
    Grad_Post_AnDA_nadirswot    = Gradient(Post_AnDA_nadirswot,2)
    FP_ConvAE_nadirswot         = itrp_FP_ConvAE_nadirswot[i,:indLat,:indLon]
    Grad_FP_ConvAE_nadirswot    = Gradient(FP_ConvAE_nadirswot,2)
    FP_GENN_nadirswot           = itrp_FP_GENN_nadirswot[i,:indLat,:indLon]
    Grad_FP_GENN_nadirswot      = Gradient(FP_GENN_nadirswot,2)
    rFP_ConvAE_nadirswot        = rec_FP_ConvAE_nadirswot[i,:indLat,:indLon]
    rGrad_FP_ConvAE_nadirswot   = Gradient(rFP_ConvAE_nadirswot,2)
    rFP_GENN_nadirswot          = rec_FP_GENN_nadirswot[i,:indLat,:indLon]
    rGrad_FP_GENN_nadirswot     = Gradient(rFP_GENN_nadirswot,2)

    ## Compute spatial coverage
    nadcov[i]		= len(np.argwhere(np.isfinite(obs_nadir.flatten())))/len(obs_nadir.flatten())
    nadswotcov[i]	= len(np.argwhere(np.isfinite(obs_nadirswot.flatten())))/len(obs_nadirswot.flatten())

    ## Compute NRMSE statistics (i.e. RMSE/stdev(gt))
    nrmse_OI_nadir[i]		= (np.sqrt(np.nanmean(((gt-np.nanmean(gt))-(OI_nadir-np.nanmean(OI_nadir)))**2)))/np.nanstd(gt)
    nrmse_Post_AnDA_nadir[i]	= (np.sqrt(np.nanmean(((gt-np.nanmean(gt))-(Post_AnDA_nadir-np.nanmean(Post_AnDA_nadir)))**2)))/np.nanstd(gt)
    nrmse_VE_DINEOF_nadir[i]	= (np.sqrt(np.nanmean(((gt-np.nanmean(gt))-(VE_DINEOF_nadir-np.nanmean(VE_DINEOF_nadir)))**2)))/np.nanstd(gt)
    nrmse_FP_ConvAE_nadir[i]    = (np.sqrt(np.nanmean(((gt-np.nanmean(gt))-(FP_ConvAE_nadir-np.nanmean(FP_ConvAE_nadir)))**2)))/np.nanstd(gt)
    nrmse_FP_GENN_nadir[i]      = (np.sqrt(np.nanmean(((gt-np.nanmean(gt))-(FP_GENN_nadir-np.nanmean(FP_GENN_nadir)))**2)))/np.nanstd(gt)
    nrmse_OI_nadirswot[i]           = (np.sqrt(np.nanmean(((gt-np.nanmean(gt))-(OI_nadirswot-np.nanmean(OI_nadirswot)))**2)))/np.nanstd(gt)
    nrmse_Post_AnDA_nadirswot[i]= (np.sqrt(np.nanmean(((gt-np.nanmean(gt))-(Post_AnDA_nadirswot-np.nanmean(Post_AnDA_nadirswot)))**2)))/np.nanstd(gt)
    nrmse_VE_DINEOF_nadirswot[i]= (np.sqrt(np.nanmean(((gt-np.nanmean(gt))-(VE_DINEOF_nadirswot-np.nanmean(VE_DINEOF_nadirswot)))**2)))/np.nanstd(gt) 
    nrmse_FP_ConvAE_nadirswot[i]  = (np.sqrt(np.nanmean(((gt-np.nanmean(gt))-(FP_ConvAE_nadirswot-np.nanmean(FP_ConvAE_nadirswot)))**2)))/np.nanstd(gt)
    nrmse_FP_GENN_nadirswot[i]    = (np.sqrt(np.nanmean(((gt-np.nanmean(gt))-(FP_GENN_nadirswot-np.nanmean(FP_GENN_nadirswot)))**2)))/np.nanstd(gt)
    ## Compute R scores 
    R_OI_nadir[i]           = 100*(1-np.nanmean(((mask1_nadir*gt-np.nanmean(mask1_nadir*gt))-(mask1_nadir*OI_nadir-np.nanmean(mask1_nadir*OI_nadir)))**2)/np.nanvar(mask1_nadir*gt))
    R_Post_AnDA_nadir[i]    = 100*(1-np.nanmean(((mask1_nadir*gt-np.nanmean(mask1_nadir*gt))-(mask1_nadir*Post_AnDA_nadir-np.nanmean(mask1_nadir*Post_AnDA_nadir)))**2)/np.nanvar(mask1_nadir*gt))
    R_VE_DINEOF_nadir[i]    = 100*(1-np.nanmean(((mask1_nadir*gt-np.nanmean(mask1_nadir*gt))-(mask1_nadir*VE_DINEOF_nadir-np.nanmean(mask1_nadir*VE_DINEOF_nadir)))**2)/np.nanvar(mask1_nadir*gt))
    R_FP_ConvAE_nadir[i]    = 100*(1-np.nanmean(((mask1_nadir*gt-np.nanmean(mask1_nadir*gt))-(mask1_nadir*FP_ConvAE_nadir-np.nanmean(mask1_nadir*FP_ConvAE_nadir)))**2)/np.nanvar(mask1_nadir*gt))
    R_FP_GENN_nadir[i]      = 100*(1-np.nanmean(((mask1_nadir*gt-np.nanmean(mask1_nadir*gt))-(mask1_nadir*FP_GENN_nadir-np.nanmean(mask1_nadir*FP_GENN_nadir)))**2)/np.nanvar(mask1_nadir*gt))
    R_OI_nadirswot[i]        = 100*(1-np.nanmean(((mask1_nadirswot*gt-np.nanmean(mask1_nadir*gt))-(mask1_nadirswot*OI_nadir-np.nanmean(mask1_nadirswot*OI_nadir)))**2)/np.nanvar(mask1_nadirswot*gt))
    R_Post_AnDA_nadirswot[i]= 100*(1-np.nanmean(((mask1_nadirswot*gt-np.nanmean(mask1_nadirswot*gt))-(mask1_nadirswot*Post_AnDA_nadirswot-np.nanmean(mask1_nadirswot*Post_AnDA_nadirswot)))**2)/np.nanvar(mask1_nadirswot*gt))
    R_VE_DINEOF_nadirswot[i]= 100*(1-np.nanmean(((mask1_nadirswot*gt-np.nanmean(mask1_nadirswot*gt))-(mask1_nadirswot*VE_DINEOF_nadirswot-np.nanmean(mask1_nadirswot*VE_DINEOF_nadirswot)))**2)/np.nanvar(mask1_nadirswot*gt))
    R_FP_ConvAE_nadirswot[i]  = 100*(1-np.nanmean(((mask1_nadirswot*gt-np.nanmean(mask1_nadirswot*gt))-(mask1_nadirswot*FP_ConvAE_nadirswot-np.nanmean(mask1_nadirswot*FP_ConvAE_nadirswot)))**2)/np.nanvar(mask1_nadirswot*gt))
    R_FP_GENN_nadirswot[i]    = 100*(1-np.nanmean(((mask1_nadirswot*gt-np.nanmean(mask1_nadirswot*gt))-(mask1_nadirswot*FP_GENN_nadirswot-np.nanmean(mask1_nadirswot*FP_GENN_nadirswot)))**2)/np.nanvar(mask1_nadirswot*gt))
    ## Compute I scores 
    I_OI_nadir[i]           = 100*(1-np.nanmean(((mask2_nadir*gt-np.nanmean(mask2_nadir*gt))-(mask2_nadir*OI_nadir-np.nanmean(mask2_nadir*OI_nadir)))**2)/np.nanvar(mask2_nadir*gt))
    I_Post_AnDA_nadir[i]    = 100*(1-np.nanmean(((mask2_nadir*gt-np.nanmean(mask2_nadir*gt))-(mask2_nadir*Post_AnDA_nadir-np.nanmean(mask2_nadir*Post_AnDA_nadir)))**2)/np.nanvar(mask2_nadir*gt))
    I_VE_DINEOF_nadir[i]    = 100*(1-np.nanmean(((mask2_nadir*gt-np.nanmean(mask2_nadir*gt))-(mask2_nadir*VE_DINEOF_nadir-np.nanmean(mask2_nadir*VE_DINEOF_nadir)))**2)/np.nanvar(mask2_nadir*gt))
    I_FP_ConvAE_nadir[i]    = 100*(1-np.nanmean(((mask2_nadir*gt-np.nanmean(mask2_nadir*gt))-(mask2_nadir*FP_ConvAE_nadir-np.nanmean(mask2_nadir*FP_ConvAE_nadir)))**2)/np.nanvar(mask2_nadir*gt))
    I_FP_GENN_nadir[i]      = 100*(1-np.nanmean(((mask2_nadir*gt-np.nanmean(mask2_nadir*gt))-(mask2_nadir*FP_GENN_nadir-np.nanmean(mask2_nadir*FP_GENN_nadir)))**2)/np.nanvar(mask2_nadir*gt))
    I_OI_nadirswot[i]        = 100*(1-np.nanmean(((mask2_nadirswot*gt-np.nanmean(mask2_nadir*gt))-(mask2_nadirswot*OI_nadir-np.nanmean(mask2_nadirswot*OI_nadir)))**2)/np.nanvar(mask2_nadirswot*gt))
    I_Post_AnDA_nadirswot[i]= 100*(1-np.nanmean(((mask2_nadirswot*gt-np.nanmean(mask2_nadirswot*gt))-(mask2_nadirswot*Post_AnDA_nadirswot-np.nanmean(mask2_nadirswot*Post_AnDA_nadirswot)))**2)/np.nanvar(mask2_nadirswot*gt))
    I_VE_DINEOF_nadirswot[i]= 100*(1-np.nanmean(((mask2_nadirswot*gt-np.nanmean(mask2_nadirswot*gt))-(mask2_nadirswot*VE_DINEOF_nadirswot-np.nanmean(mask2_nadirswot*VE_DINEOF_nadirswot)))**2)/np.nanvar(mask2_nadirswot*gt))
    I_FP_ConvAE_nadirswot[i]  = 100*(1-np.nanmean(((mask2_nadirswot*gt-np.nanmean(mask2_nadirswot*gt))-(mask2_nadirswot*FP_ConvAE_nadirswot-np.nanmean(mask2_nadirswot*FP_ConvAE_nadirswot)))**2)/np.nanvar(mask2_nadirswot*gt))
    I_FP_GENN_nadirswot[i]    = 100*(1-np.nanmean(((mask2_nadirswot*gt-np.nanmean(mask2_nadirswot*gt))-(mask2_nadirswot*FP_GENN_nadirswot-np.nanmean(mask2_nadirswot*FP_GENN_nadirswot)))**2)/np.nanvar(mask2_nadirswot*gt))
    ## Compute AE scores 
    AE_FP_ConvAE_nadir[i]    = 100*(1-np.nanmean(((gt-np.nanmean(gt))-(rFP_ConvAE_nadir-np.nanmean(rFP_ConvAE_nadir)))**2)/np.nanvar(gt))
    AE_FP_GENN_nadir[i]      = 100*(1-np.nanmean(((gt-np.nanmean(gt))-(rFP_GENN_nadir-np.nanmean(rFP_GENN_nadir)))**2)/np.nanvar(gt))
    AE_FP_ConvAE_nadirswot[i]  = 100*(1-np.nanmean(((gt-np.nanmean(gt))-(rFP_ConvAE_nadirswot-np.nanmean(rFP_ConvAE_nadirswot)))**2)/np.nanvar(gt))
    AE_FP_GENN_nadirswot[i]    = 100*(1-np.nanmean(((gt-np.nanmean(gt))-(rFP_GENN_nadirswot-np.nanmean(rFP_GENN_nadirswot)))**2)/np.nanvar(gt))

    ## Compute NRMSE Grad statistics 
    nrmse_Grad_OI_nadir[i]		= (np.sqrt(np.nanmean(((Grad_gt-np.nanmean(Grad_gt))-(Grad_OI_nadir-np.nanmean(Grad_OI_nadir)))**2)))/np.nanstd(Grad_gt)
    nrmse_Grad_Post_AnDA_nadir[i]	= (np.sqrt(np.nanmean(((Grad_gt-np.nanmean(Grad_gt))-(Grad_Post_AnDA_nadir-np.nanmean(Grad_Post_AnDA_nadir)))**2)))/np.nanstd(Grad_gt)
    nrmse_Grad_VE_DINEOF_nadir[i]	= (np.sqrt(np.nanmean(((Grad_gt-np.nanmean(Grad_gt))-(Grad_VE_DINEOF_nadir-np.nanmean(Grad_VE_DINEOF_nadir)))**2)))/np.nanstd(Grad_gt)
    nrmse_Grad_FP_ConvAE_nadir[i]    = (np.sqrt(np.nanmean(((Grad_gt-np.nanmean(Grad_gt))-(Grad_FP_ConvAE_nadir-np.nanmean(Grad_FP_ConvAE_nadir)))**2)))/np.nanstd(Grad_gt)
    nrmse_Grad_FP_GENN_nadir[i]      = (np.sqrt(np.nanmean(((Grad_gt-np.nanmean(Grad_gt))-(Grad_FP_GENN_nadir-np.nanmean(Grad_FP_GENN_nadir)))**2)))/np.nanstd(Grad_gt)
    nrmse_Grad_OI_nadirswot[i]           = (np.sqrt(np.nanmean(((Grad_gt-np.nanmean(Grad_gt))-(Grad_OI_nadirswot-np.nanmean(Grad_OI_nadirswot)))**2)))/np.nanstd(Grad_gt)
    nrmse_Grad_Post_AnDA_nadirswot[i]= (np.sqrt(np.nanmean(((Grad_gt-np.nanmean(Grad_gt))-(Grad_Post_AnDA_nadirswot-np.nanmean(Grad_Post_AnDA_nadirswot)))**2)))/np.nanstd(Grad_gt)
    nrmse_Grad_VE_DINEOF_nadirswot[i]= (np.sqrt(np.nanmean(((Grad_gt-np.nanmean(Grad_gt))-(Grad_VE_DINEOF_nadirswot-np.nanmean(Grad_VE_DINEOF_nadirswot)))**2)))/np.nanstd(Grad_gt) 
    nrmse_Grad_FP_ConvAE_nadirswot[i]  = (np.sqrt(np.nanmean(((Grad_gt-np.nanmean(Grad_gt))-(Grad_FP_ConvAE_nadirswot-np.nanmean(Grad_FP_ConvAE_nadirswot)))**2)))/np.nanstd(Grad_gt)
    nrmse_Grad_FP_GENN_nadirswot[i]    = (np.sqrt(np.nanmean(((Grad_gt-np.nanmean(Grad_gt))-(Grad_FP_GENN_nadirswot-np.nanmean(Grad_FP_GENN_nadirswot)))**2)))/np.nanstd(Grad_gt)
    ## Compute R scores 
    R_Grad_OI_nadir[i]           = 100*(1-np.nanmean(((mask1_nadir*Grad_gt-np.nanmean(mask1_nadir*Grad_gt))-(mask1_nadir*Grad_OI_nadir-np.nanmean(mask1_nadir*Grad_OI_nadir)))**2)/np.nanvar(mask1_nadir*Grad_gt))
    R_Grad_Post_AnDA_nadir[i]    = 100*(1-np.nanmean(((mask1_nadir*Grad_gt-np.nanmean(mask1_nadir*Grad_gt))-(mask1_nadir*Grad_Post_AnDA_nadir-np.nanmean(mask1_nadir*Grad_Post_AnDA_nadir)))**2)/np.nanvar(mask1_nadir*Grad_gt))
    R_Grad_VE_DINEOF_nadir[i]    = 100*(1-np.nanmean(((mask1_nadir*Grad_gt-np.nanmean(mask1_nadir*Grad_gt))-(mask1_nadir*Grad_VE_DINEOF_nadir-np.nanmean(mask1_nadir*Grad_VE_DINEOF_nadir)))**2)/np.nanvar(mask1_nadir*Grad_gt))
    R_Grad_FP_ConvAE_nadir[i]    = 100*(1-np.nanmean(((mask1_nadir*Grad_gt-np.nanmean(mask1_nadir*Grad_gt))-(mask1_nadir*Grad_FP_ConvAE_nadir-np.nanmean(mask1_nadir*Grad_FP_ConvAE_nadir)))**2)/np.nanvar(mask1_nadir*Grad_gt))
    R_Grad_FP_GENN_nadir[i]      = 100*(1-np.nanmean(((mask1_nadir*Grad_gt-np.nanmean(mask1_nadir*Grad_gt))-(mask1_nadir*Grad_FP_GENN_nadir-np.nanmean(mask1_nadir*Grad_FP_GENN_nadir)))**2)/np.nanvar(mask1_nadir*Grad_gt))
    R_Grad_OI_nadirswot[i]        = 100*(1-np.nanmean(((mask1_nadirswot*Grad_gt-np.nanmean(mask1_nadir*Grad_gt))-(mask1_nadirswot*Grad_OI_nadir-np.nanmean(mask1_nadirswot*Grad_OI_nadir)))**2)/np.nanvar(mask1_nadirswot*Grad_gt))
    R_Grad_Post_AnDA_nadirswot[i]= 100*(1-np.nanmean(((mask1_nadirswot*Grad_gt-np.nanmean(mask1_nadirswot*Grad_gt))-(mask1_nadirswot*Grad_Post_AnDA_nadirswot-np.nanmean(mask1_nadirswot*Grad_Post_AnDA_nadirswot)))**2)/np.nanvar(mask1_nadirswot*Grad_gt))
    R_Grad_VE_DINEOF_nadirswot[i]= 100*(1-np.nanmean(((mask1_nadirswot*Grad_gt-np.nanmean(mask1_nadirswot*Grad_gt))-(mask1_nadirswot*Grad_VE_DINEOF_nadirswot-np.nanmean(mask1_nadirswot*Grad_VE_DINEOF_nadirswot)))**2)/np.nanvar(mask1_nadirswot*Grad_gt))
    R_Grad_FP_ConvAE_nadirswot[i]  = 100*(1-np.nanmean(((mask1_nadirswot*Grad_gt-np.nanmean(mask1_nadirswot*Grad_gt))-(mask1_nadirswot*Grad_FP_ConvAE_nadirswot-np.nanmean(mask1_nadirswot*Grad_FP_ConvAE_nadirswot)))**2)/np.nanvar(mask1_nadirswot*Grad_gt))
    R_Grad_FP_GENN_nadirswot[i]    = 100*(1-np.nanmean(((mask1_nadirswot*Grad_gt-np.nanmean(mask1_nadirswot*Grad_gt))-(mask1_nadirswot*Grad_FP_GENN_nadirswot-np.nanmean(mask1_nadirswot*Grad_FP_GENN_nadirswot)))**2)/np.nanvar(mask1_nadirswot*Grad_gt))
    ## Compute I scores 
    I_Grad_OI_nadir[i]           = 100*(1-np.nanmean(((mask2_nadir*Grad_gt-np.nanmean(mask2_nadir*Grad_gt))-(mask2_nadir*Grad_OI_nadir-np.nanmean(mask2_nadir*Grad_OI_nadir)))**2)/np.nanvar(mask2_nadir*Grad_gt))
    I_Grad_Post_AnDA_nadir[i]    = 100*(1-np.nanmean(((mask2_nadir*Grad_gt-np.nanmean(mask2_nadir*Grad_gt))-(mask2_nadir*Grad_Post_AnDA_nadir-np.nanmean(mask2_nadir*Grad_Post_AnDA_nadir)))**2)/np.nanvar(mask2_nadir*Grad_gt))
    I_Grad_VE_DINEOF_nadir[i]    = 100*(1-np.nanmean(((mask2_nadir*Grad_gt-np.nanmean(mask2_nadir*Grad_gt))-(mask2_nadir*Grad_VE_DINEOF_nadir-np.nanmean(mask2_nadir*Grad_VE_DINEOF_nadir)))**2)/np.nanvar(mask2_nadir*Grad_gt))
    I_Grad_FP_ConvAE_nadir[i]    = 100*(1-np.nanmean(((mask2_nadir*Grad_gt-np.nanmean(mask2_nadir*Grad_gt))-(mask2_nadir*Grad_FP_ConvAE_nadir-np.nanmean(mask2_nadir*Grad_FP_ConvAE_nadir)))**2)/np.nanvar(mask2_nadir*Grad_gt))
    I_Grad_FP_GENN_nadir[i]      = 100*(1-np.nanmean(((mask2_nadir*Grad_gt-np.nanmean(mask2_nadir*Grad_gt))-(mask2_nadir*Grad_FP_GENN_nadir-np.nanmean(mask2_nadir*Grad_FP_GENN_nadir)))**2)/np.nanvar(mask2_nadir*Grad_gt))
    I_Grad_OI_nadirswot[i]        = 100*(1-np.nanmean(((mask2_nadirswot*Grad_gt-np.nanmean(mask2_nadir*Grad_gt))-(mask2_nadirswot*Grad_OI_nadir-np.nanmean(mask2_nadirswot*Grad_OI_nadir)))**2)/np.nanvar(mask2_nadirswot*Grad_gt))
    I_Grad_Post_AnDA_nadirswot[i]= 100*(1-np.nanmean(((mask2_nadirswot*Grad_gt-np.nanmean(mask2_nadirswot*Grad_gt))-(mask2_nadirswot*Grad_Post_AnDA_nadirswot-np.nanmean(mask2_nadirswot*Grad_Post_AnDA_nadirswot)))**2)/np.nanvar(mask2_nadirswot*Grad_gt))
    I_Grad_VE_DINEOF_nadirswot[i]= 100*(1-np.nanmean(((mask2_nadirswot*Grad_gt-np.nanmean(mask2_nadirswot*Grad_gt))-(mask2_nadirswot*Grad_VE_DINEOF_nadirswot-np.nanmean(mask2_nadirswot*Grad_VE_DINEOF_nadirswot)))**2)/np.nanvar(mask2_nadirswot*Grad_gt))
    I_Grad_FP_ConvAE_nadirswot[i]  = 100*(1-np.nanmean(((mask2_nadirswot*Grad_gt-np.nanmean(mask2_nadirswot*Grad_gt))-(mask2_nadirswot*Grad_FP_ConvAE_nadirswot-np.nanmean(mask2_nadirswot*Grad_FP_ConvAE_nadirswot)))**2)/np.nanvar(mask2_nadirswot*Grad_gt))
    I_Grad_FP_GENN_nadirswot[i]    = 100*(1-np.nanmean(((mask2_nadirswot*Grad_gt-np.nanmean(mask2_nadirswot*Grad_gt))-(mask2_nadirswot*Grad_FP_GENN_nadirswot-np.nanmean(mask2_nadirswot*Grad_FP_GENN_nadirswot)))**2)/np.nanvar(mask2_nadirswot*Grad_gt))
    ## Compute AE scores 
    AE_Grad_FP_ConvAE_nadir[i]    = 100*(1-np.nanmean(((Grad_gt-np.nanmean(Grad_gt))-(rGrad_FP_ConvAE_nadir-np.nanmean(rGrad_FP_ConvAE_nadir)))**2)/np.nanvar(Grad_gt))
    AE_Grad_FP_GENN_nadir[i]      = 100*(1-np.nanmean(((Grad_gt-np.nanmean(Grad_gt))-(rGrad_FP_GENN_nadir-np.nanmean(rGrad_FP_GENN_nadir)))**2)/np.nanvar(Grad_gt))
    AE_Grad_FP_ConvAE_nadirswot[i]  = 100*(1-np.nanmean(((Grad_gt-np.nanmean(Grad_gt))-(rGrad_FP_ConvAE_nadirswot-np.nanmean(rGrad_FP_ConvAE_nadirswot)))**2)/np.nanvar(Grad_gt))
    AE_Grad_FP_GENN_nadirswot[i]    = 100*(1-np.nanmean(((Grad_gt-np.nanmean(Grad_gt))-(rGrad_FP_GENN_nadirswot-np.nanmean(rGrad_FP_GENN_nadirswot)))**2)/np.nanvar(Grad_gt))

    '''# Display individual maps
    var=['gt','obs_nadir','obs_nadirswot',\
         'OI_nadir','OI_nadirswot',\
         'AnDA_nadir','AnDA_nadirswot',\
         'Post_AnDA_nadir','Post_AnDA_nadirswot',\
         'FP_ConvAE_nadir','FP_ConvAE_nadirswot',\
         'FP_GENN_nadir','FP_GENN_nadirswot',\
         'VE_DINEOF_nadir','VE_DINEOF_nadirswot',]
    title=['Ground Truth','Obs (nadir)','Obs (nadir+swot)',\
           'OI (nadir)','OI (nadir+swot)',\
           'AnDA (nadir)','AnDA (nadir+swot)',\
           'Post-AnDA (nadir)','Post-AnDA (nadir+swot)',\
           'FP-ConvAE (nadir)','FP-ConvAE (nadir+swot)',\
           'FP-GENN (nadir)','FP-GENN (nadir+swot)',\
           'VE-DINEOF (nadir)','VE-DINEOF (nadir+swot)']
    for ivar in range(len(var)):
        resfile = workpath+"/results_"+var[ivar]+'_'+day+".png"
        fig, ax = plt.subplots(1,1,figsize=(10,10),
                          subplot_kw=dict(projection=ccrs.PlateCarree(central_longitude=0.0)))   
        vmax = np.nanmax(np.abs(gt))
        vmin = -1.*vmax
        cmap="coolwarm"
        plot2(ax,lon,lat,eval(var[ivar]),title[ivar],\
             extent=extent,cmap=cmap,vmin=vmin,vmax=vmax)
        plt.savefig(resfile)       # save the figure
        plt.close()                 # close the figure

    # Display individual gradient maps
    var=['Grad_gt',\
         'Grad_OI_nadir','Grad_OI_nadirswot',\
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
    for ivar in range(len(var)):
        resfile = workpath+"/results_"+var[ivar]+'_'+day+".png"
        fig, ax = plt.subplots(1,1,figsize=(10,10),
                          subplot_kw=dict(projection=ccrs.PlateCarree(central_longitude=0.0)))
        vmin = 0 ; vmax = np.nanmax(Grad_gt)
        cmap="viridis"
        plot2(ax,lon,lat,eval(var[ivar]),title[ivar],\
             extent=extent,cmap=cmap,vmin=vmin,vmax=vmax)
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
             extent=extent,cmap=cmap,vmin=vmin,vmax=vmax)
    ax[1][0].set_visible(False)
    for ivar in range(0,len(var)):
        i = int(np.floor(ivar/2))+1 ; j = (ivar%2)
        plot(ax,j,i,lon,lat,eval(var[ivar]),title[ivar],\
             extent=extent,cmap=cmap,vmin=vmin,vmax=vmax)
    plt.subplots_adjust(hspace=0.3,wspace=0.5)
    plt.savefig(resfile1)	# save the figure
    plt.close()			# close the figure

    # Display gradients
    var=['obs_nadir','obs_nadirswot',\
         'Grad_OI_nadir','Grad_OI_nadirswot',\
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
             extent=extent,cmap=cmap,vmin=vmin,vmax=vmax)
    ax[1][0].set_visible(False)
    for ivar in range(0,len(var)):
        i = int(np.floor(ivar/2))+1 ; j = (ivar%2)
        if i==1:
            plot(ax,j,i,lon,lat,eval(var[ivar]),title[ivar],\
                 extent=extent,cmap=cmap2,vmin=vmin2,vmax=vmax2)
        else:
            plot(ax,j,i,lon,lat,eval(var[ivar]),title[ivar],\
                 extent=extent,cmap=cmap,vmin=vmin,vmax=vmax)
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
                styles=['k','p','o','p','o','p','o','p','o','p','o','p','o'],\
                colors=['k','y','y','mediumseagreen','mediumseagreen',\
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
                styles=['k','p','o','p','o','p','o','p','o','p','o','p','o'],\
                colors=['k','y','y','mediumseagreen','mediumseagreen',\
                                'seagreen','seagreen',\
                                'steelblue','steelblue',\
                                'mediumorchid','mediumorchid',\
                                'darkorange','darkorange'])
    plt.savefig(resfile)
    plt.close()

    ## Radial Power Spectrum (RAPS)
    resfile=workpath+"/results_RAPS_"+day+".png"
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
    resfile=workpath+"/results_diff_RAPS_"+day+".png"
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
    resfile=workpath+"/results_diff_avg_RAPS_nadir"+day+".png"
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
    resfile=workpath+"/results_diff_avg_RAPS_nadirswot"+day+".png"
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

## SSH score tables (mean of daily scores)
index=list(range(5,16))
index.extend(range(25,36))
index.extend(range(45,56))
index.extend(range(65,76))
tab_scores = np.zeros((10,3))
tab_scores[0,0] = np.nanmean(R_OI_nadir)
tab_scores[0,1] = np.nanmean(I_OI_nadir)
tab_scores[0,2] = np.nan
tab_scores[1,0] = np.nanmean(R_Post_AnDA_nadir)
tab_scores[1,1] = np.nanmean(I_Post_AnDA_nadir)
tab_scores[1,2] = np.nan
tab_scores[2,0] = np.nanmean(R_VE_DINEOF_nadir)
tab_scores[2,1] = np.nanmean(I_VE_DINEOF_nadir)
tab_scores[2,2] = np.nan
tab_scores[3,0] = np.nanmean(R_FP_ConvAE_nadir)
tab_scores[3,1] = np.nanmean(I_FP_ConvAE_nadir)
tab_scores[3,2] = np.nanmean(AE_FP_ConvAE_nadir)
tab_scores[4,0] = np.nanmean(R_FP_GENN_nadir)
tab_scores[4,1] = np.nanmean(I_FP_GENN_nadir)
tab_scores[4,2] = np.nanmean(AE_FP_GENN_nadir)
tab_scores[5,0] = np.nanmean(R_OI_nadirswot)
tab_scores[5,1] = np.nanmean(I_OI_nadirswot)
tab_scores[5,2] = np.nan
tab_scores[6,0] = np.nanmean(R_Post_AnDA_nadirswot)
tab_scores[6,1] = np.nanmean(I_Post_AnDA_nadirswot)
tab_scores[6,2] = np.nan
tab_scores[7,0] = np.nanmean(R_VE_DINEOF_nadirswot)
tab_scores[7,1] = np.nanmean(I_VE_DINEOF_nadirswot)
tab_scores[7,2] = np.nan
tab_scores[8,0] = np.nanmean(R_FP_ConvAE_nadirswot)
tab_scores[8,1] = np.nanmean(I_FP_ConvAE_nadirswot)
tab_scores[8,2] = np.nanmean(AE_FP_ConvAE_nadirswot)
tab_scores[9,0] = np.nanmean(R_FP_GENN_nadirswot)
tab_scores[9,1] = np.nanmean(I_FP_GENN_nadirswot)
tab_scores[9,2] = np.nanmean(AE_FP_GENN_nadirswot)
np.savetxt(fname=workpath+"/tab_scores_SSH.txt",X=tab_scores,fmt='%2.2f')

## GradSSH score tables (mean of daily scores)
tab_scores = np.zeros((10,3))
tab_scores[0,0] = np.nanmean(R_Grad_OI_nadir)
tab_scores[0,1] = np.nanmean(I_Grad_OI_nadir)
tab_scores[0,2] = np.nan
tab_scores[1,0] = np.nanmean(R_Grad_Post_AnDA_nadir)
tab_scores[1,1] = np.nanmean(I_Grad_Post_AnDA_nadir)
tab_scores[1,2] = np.nan
tab_scores[2,0] = np.nanmean(R_Grad_VE_DINEOF_nadir)
tab_scores[2,1] = np.nanmean(I_Grad_VE_DINEOF_nadir)
tab_scores[2,2] = np.nan
tab_scores[3,0] = np.nanmean(R_Grad_FP_ConvAE_nadir)
tab_scores[3,1] = np.nanmean(I_Grad_FP_ConvAE_nadir)
tab_scores[3,2] = np.nanmean(AE_Grad_FP_ConvAE_nadir)
tab_scores[4,0] = np.nanmean(R_Grad_FP_GENN_nadir)
tab_scores[4,1] = np.nanmean(I_Grad_FP_GENN_nadir)
tab_scores[4,2] = np.nanmean(AE_Grad_FP_GENN_nadir)
tab_scores[5,0] = np.nanmean(R_Grad_OI_nadirswot)
tab_scores[5,1] = np.nanmean(I_Grad_OI_nadirswot)
tab_scores[5,2] = np.nan
tab_scores[6,0] = np.nanmean(R_Grad_Post_AnDA_nadirswot)
tab_scores[6,1] = np.nanmean(I_Grad_Post_AnDA_nadirswot)
tab_scores[6,2] = np.nan
tab_scores[7,0] = np.nanmean(R_Grad_VE_DINEOF_nadirswot)
tab_scores[7,1] = np.nanmean(I_Grad_VE_DINEOF_nadirswot)
tab_scores[7,2] = np.nan
tab_scores[8,0] = np.nanmean(R_Grad_FP_ConvAE_nadirswot)
tab_scores[8,1] = np.nanmean(I_Grad_FP_ConvAE_nadirswot)
tab_scores[8,2] = np.nanmean(AE_Grad_FP_ConvAE_nadirswot)
tab_scores[9,0] = np.nanmean(R_Grad_FP_GENN_nadirswot)
tab_scores[9,1] = np.nanmean(I_Grad_FP_GENN_nadirswot)
tab_scores[9,2] = np.nanmean(AE_Grad_FP_GENN_nadirswot)
np.savetxt(fname=workpath+"/tab_scores_GradSSH.txt",X=tab_scores,fmt='%2.2f')

def Iscore(mask1,gt,itrp):
    return 100*(1-np.nanmean(((mask1*gt-np.nanmean(mask1*gt))-(mask1*itrp-np.nanmean(mask1*itrp)))**2)/np.nanvar(mask1*gt))
def Rscore(mask1,gt,itrp):
    return 100*(1-np.nanmean(((mask1*gt-np.nanmean(mask1*gt))-(mask1*itrp-np.nanmean(mask1*itrp)))**2)/np.nanvar(mask1*gt))
def AEscore(gt,itrp):
    return 100*(1-np.nanmean(((gt-np.nanmean(gt))-(itrp-np.nanmean(itrp)))**2)/np.nanvar(gt))
obs_nadir   = AnDA_ssh_1_nadir.Obs[:,:indLat,:indLon].flatten()
mask1_nadir = np.where(np.isnan(obs_nadir),np.nan,1)
mask2_nadir = np.where(np.isnan(obs_nadir),1,np.nan)
obs_nadirswot   = AnDA_ssh_1_nadirswot.Obs[:,:indLat,:indLon].flatten()
mask1_nadirswot = np.where(np.isnan(obs_nadirswot),np.nan,1)
mask2_nadirswot = np.where(np.isnan(obs_nadirswot),1,np.nan)

## SSH score tables
index=list(range(5,16))
index.extend(range(25,36))
index.extend(range(45,56))
index.extend(range(65,76))
# apply HPF to visualize Taylor diagrams only for small scales
HR = AnDA_ssh_1.itrp_OI[:,:indLat,:indLon]
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
# create scores table
tab_scores = np.zeros((10,3))
tab_scores[0,0] = Rscore(mask1_nadir,AnDA_ssh_1.GT[:,:indLat,:indLon].flatten()-lr,AnDA_ssh_1_nadir.itrp_OI[:,:indLat,:indLon].flatten()-lr)
tab_scores[0,1] = Iscore(mask2_nadir,AnDA_ssh_1.GT[:,:indLat,:indLon].flatten()-lr,AnDA_ssh_1_nadir.itrp_OI[:,:indLat,:indLon].flatten()-lr)
tab_scores[0,2] = np.nan
tab_scores[1,0] = Rscore(mask1_nadir,AnDA_ssh_1.GT[:,:indLat,:indLon].flatten()-lr,AnDA_ssh_1_nadir.itrp_postAnDA[:,:indLat,:indLon].flatten()-lr)
tab_scores[1,1] = Iscore(mask2_nadir,AnDA_ssh_1.GT[:,:indLat,:indLon].flatten()-lr,AnDA_ssh_1_nadir.itrp_postAnDA[:,:indLat,:indLon].flatten()-lr)
tab_scores[1,2] = np.nan
tab_scores[2,0] = Rscore(mask1_nadir,AnDA_ssh_1.GT[:,:indLat,:indLon].flatten()-lr,itrp_dineof_nadir[:,:indLat,:indLon].flatten()-lr)
tab_scores[2,1] = Iscore(mask2_nadir,AnDA_ssh_1.GT[:,:indLat,:indLon].flatten()-lr,itrp_dineof_nadir[:,:indLat,:indLon].flatten()-lr)
tab_scores[2,2] = np.nan
tab_scores[3,0] = Rscore(mask1_nadir,AnDA_ssh_1.GT[:,:indLat,:indLon].flatten()-lr,itrp_FP_ConvAE_nadir[:,:indLat,:indLon].flatten()-lr)
tab_scores[3,1] = Iscore(mask2_nadir,AnDA_ssh_1.GT[:,:indLat,:indLon].flatten()-lr,itrp_FP_ConvAE_nadir[:,:indLat,:indLon].flatten()-lr)
tab_scores[3,2] = AEscore(AnDA_ssh_1.GT[:,:indLat,:indLon].flatten()-lr,rec_FP_ConvAE_nadir[:,:indLat,:indLon].flatten()-lr)
tab_scores[4,0] = Rscore(mask1_nadir,AnDA_ssh_1.GT[:,:indLat,:indLon].flatten()-lr,itrp_FP_GENN_nadir[:,:indLat,:indLon].flatten()-lr)
tab_scores[4,1] = Iscore(mask2_nadir,AnDA_ssh_1.GT[:,:indLat,:indLon].flatten()-lr,itrp_FP_GENN_nadir[:,:indLat,:indLon].flatten()-lr)
tab_scores[4,2] = AEscore(AnDA_ssh_1.GT[:,:indLat,:indLon].flatten()-lr,rec_FP_GENN_nadir[:,:indLat,:indLon].flatten()-lr)
tab_scores[5,0] = Rscore(mask1_nadirswot,AnDA_ssh_1.GT[:,:indLat,:indLon].flatten()-lr,AnDA_ssh_1_nadirswot.itrp_OI[:,:indLat,:indLon].flatten()-lr)
tab_scores[5,1] = Iscore(mask2_nadirswot,AnDA_ssh_1.GT[:,:indLat,:indLon].flatten()-lr,AnDA_ssh_1_nadirswot.itrp_OI[:,:indLat,:indLon].flatten()-lr)
tab_scores[5,2] = np.nan
tab_scores[6,0] = Rscore(mask1_nadirswot,AnDA_ssh_1.GT[:,:indLat,:indLon].flatten()-lr,AnDA_ssh_1_nadirswot.itrp_postAnDA[:,:indLat,:indLon].flatten()-lr)
tab_scores[6,1] = Iscore(mask2_nadirswot,AnDA_ssh_1.GT[:,:indLat,:indLon].flatten()-lr,AnDA_ssh_1_nadirswot.itrp_postAnDA[:,:indLat,:indLon].flatten()-lr)
tab_scores[6,2] = np.nan
tab_scores[7,0] = Rscore(mask1_nadirswot,AnDA_ssh_1.GT[:,:indLat,:indLon].flatten()-lr,itrp_dineof_nadirswot[:,:indLat,:indLon].flatten()-lr)
tab_scores[7,1] = Iscore(mask2_nadirswot,AnDA_ssh_1.GT[:,:indLat,:indLon].flatten()-lr,itrp_dineof_nadirswot[:,:indLat,:indLon].flatten()-lr)
tab_scores[7,2] = np.nan
tab_scores[8,0] = Rscore(mask1_nadirswot,AnDA_ssh_1.GT[:,:indLat,:indLon].flatten()-lr,itrp_FP_ConvAE_nadirswot[:,:indLat,:indLon].flatten()-lr)
tab_scores[8,1] = Iscore(mask2_nadirswot,AnDA_ssh_1.GT[:,:indLat,:indLon].flatten()-lr,itrp_FP_ConvAE_nadirswot[:,:indLat,:indLon].flatten()-lr)
tab_scores[8,2] = AEscore(AnDA_ssh_1.GT[:,:indLat,:indLon].flatten()-lr,rec_FP_ConvAE_nadirswot[:,:indLat,:indLon].flatten()-lr)
tab_scores[9,0] = Rscore(mask1_nadirswot,AnDA_ssh_1.GT[:,:indLat,:indLon].flatten()-lr,itrp_FP_GENN_nadirswot[:,:indLat,:indLon].flatten()-lr)
tab_scores[9,1] = Iscore(mask2_nadirswot,AnDA_ssh_1.GT[:,:indLat,:indLon].flatten()-lr,itrp_FP_GENN_nadirswot[:,:indLat,:indLon].flatten()-lr)
tab_scores[9,2] = AEscore(AnDA_ssh_1.GT[:,:indLat,:indLon].flatten()-lr,rec_FP_GENN_nadirswot[:,:indLat,:indLon].flatten()-lr)
np.savetxt(fname=workpath+"/tab_scores_SSH_2.txt",X=tab_scores,fmt='%2.2f')

## Taylor diagrams
resfile=workpath+"/Taylor_diagram.png"
label=['GT',\
       'OI (nadir)','OI (nadir+swot)',\
       'AnDA (nadir)','AnDA (nadir+swot)',\
       'Post-AnDA (nadir)','Post-AnDA (nadir+swot)',\
       'FP-ConvAE (nadir)','FP-ConvAE (nadir+swot)',\
       'FP-GENN (nadir)','FP-GENN (nadir+swot)',\
       'VE-DINEOF (nadir)','VE-DINEOF (nadir+swot)']
series={'gt':AnDA_ssh_1.GT[:,:indLat,:indLon].flatten()-lr,
        'OI_nadir':AnDA_ssh_1_nadir.itrp_OI[:,:indLat,:indLon].flatten()-lr,
        'OI_nadirswot':AnDA_ssh_1_nadirswot.itrp_OI[:,:indLat,:indLon].flatten()-lr,
        'AnDA_nadir':AnDA_ssh_1_nadir.itrp_AnDA[:,:indLat,:indLon].flatten()-lr,\
        'AnDA_nadirswot':AnDA_ssh_1_nadirswot.itrp_AnDA[:,:indLat,:indLon].flatten()-lr,
        'Post_AnDA_nadir':AnDA_ssh_1_nadir.itrp_postAnDA[:,:indLat,:indLon].flatten()-lr,\
        'Post_AnDA_nadirswot':AnDA_ssh_1_nadirswot.itrp_postAnDA[:,:indLat,:indLon].flatten()-lr,
        'FP_ConvAE_nadir':itrp_FP_ConvAE_nadir[:,:indLat,:indLon].flatten()-lr,\
        'FP_ConvAE_nadirswot':itrp_FP_ConvAE_nadirswot[:,:indLat,:indLon].flatten()-lr,
        'FP_GENN_nadir':itrp_FP_GENN_nadir[:,:indLat,:indLon].flatten()-lr,\
        'FP_GENN_nadirswot':itrp_FP_GENN_nadirswot[:,:indLat,:indLon].flatten()-lr,
        'VE_DINEOF_nadir':itrp_dineof_nadir[:,:indLat,:indLon].flatten()-lr,\
        'VE_DINEOF_nadirswot':itrp_dineof_nadirswot[:,:indLat,:indLon].flatten()-lr}
Taylor_diag(series,label,\
            styles=['k','p','o','p','o','p','o','p','o','p','o','p','o'],\
            colors=['k','y','y','mediumseagreen','mediumseagreen',\
                                'seagreen','seagreen',\
                                'steelblue','steelblue',\
                                'mediumorchid','mediumorchid',\
                                'darkorange','darkorange'])
plt.savefig(resfile)
plt.close()

## Plot averaged RAPS
resfile=workpath+"/results_avg_RAPS.png"
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
resfile=workpath+"/results_avg_RAPS_nadir.png"
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
resfile=workpath+"/results_avg_RAPS_nadirswot.png"
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
resfile=workpath+"/results_diff_avg_RAPS.png"
f0_nadir, Pf_OI_nadir                = avg_err_raPsd2dv1(AnDA_ssh_1_nadir.itrp_OI[:,:indLat,:indLon],AnDA_ssh_1.GT[:,:indLat,:indLon],resssh,True)
f0_nadirswot, Pf_OI_nadirswot        = avg_err_raPsd2dv1(AnDA_ssh_1_nadirswot.itrp_OI[:,:indLat,:indLon],AnDA_ssh_1.GT[:,:indLat,:indLon],resssh,True)
f1_nadir, Pf_AnDA_nadir              = avg_err_raPsd2dv1(AnDA_ssh_1_nadir.itrp_AnDA[:,:indLat,:indLon],AnDA_ssh_1.GT[:,:indLat,:indLon],resssh,True)
f1_nadirswot, Pf_AnDA_nadirswot      = avg_err_raPsd2dv1(AnDA_ssh_1_nadirswot.itrp_AnDA[:,:indLat,:indLon],AnDA_ssh_1.GT[:,:indLat,:indLon],resssh,True)
f2_nadir, Pf_postAnDA_nadir          = avg_err_raPsd2dv1(AnDA_ssh_1_nadir.itrp_postAnDA[:,:indLat,:indLon],AnDA_ssh_1.GT[:,:indLat,:indLon],resssh,True)
f2_nadirswot, Pf_postAnDA_nadirswot  = avg_err_raPsd2dv1(AnDA_ssh_1_nadirswot.itrp_postAnDA[:,:indLat,:indLon],AnDA_ssh_1.GT[:,:indLat,:indLon],resssh,True)
f3_nadir, Pf_VE_DINEOF_nadir         = avg_err_raPsd2dv1(itrp_dineof_nadir[:,:indLat,:indLon],AnDA_ssh_1.GT[:,:indLat,:indLon],resssh,True)
f3_nadirswot,Pf_VE_DINEOF_nadirswot  = avg_err_raPsd2dv1(itrp_dineof_nadirswot[:,:indLat,:indLon],AnDA_ssh_1.GT[:,:indLat,:indLon],resssh,True)
f4_nadir, Pf_FP_ConvAE_nadir         = avg_err_raPsd2dv1(itrp_FP_ConvAE_nadir[:,:indLat,:indLon],AnDA_ssh_1.GT[:,:indLat,:indLon],resssh,True)
f4_nadirswot, Pf_FP_ConvAE_nadirswot = avg_err_raPsd2dv1(itrp_FP_ConvAE_nadirswot[:,:indLat,:indLon],AnDA_ssh_1.GT[:,:indLat,:indLon],resssh,True)
f5_nadir, Pf_FP_GENN_nadir           = avg_err_raPsd2dv1(itrp_FP_GENN_nadir[:,:indLat,:indLon],AnDA_ssh_1.GT[:,:indLat,:indLon],resssh,True)
f5_nadirswot, Pf_FP_GENN_nadirswot   = avg_err_raPsd2dv1(itrp_FP_GENN_nadirswot[:,:indLat,:indLon],AnDA_ssh_1.GT[:,:indLat,:indLon],resssh,True)
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
resfile=workpath+"/results_diff_avg_RAPS_nadir.png"
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
resfile=workpath+"/results_diff_avg_RAPS_nadirswot.png"
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
p2 = axes2.bar(range(N), nadswotcov-nadcov, width,bottom=nadcov,color='g',alpha=0.25)
axes2.set_ylim(0, 1)
axes2.set_ylabel('Spatial Coverage (%)')
axes2.margins(x=0)
resfile=workpath+"/TS_nRMSE.png"
plt.savefig(resfile,bbox_inches="tight")    # save the figure
plt.close()         	# close the figure

## Plot time series (nadir)
N = len(lday)
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
p2 = axes2.bar(range(N), nadswotcov-nadcov, width,bottom=nadcov,color='g',alpha=0.25)
axes2.set_ylim(0, 1)
axes2.set_ylabel('Spatial Coverage (%)')
axes2.margins(x=0)
resfile=workpath+"/TS_nRMSE_nadir.png"
plt.savefig(resfile,bbox_inches="tight")    # save the figure
plt.close()             # close the figure

## Plot time series (nadir/swot)
N = len(lday)
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
p2 = axes2.bar(range(N), nadswotcov-nadcov, width,bottom=nadcov,color='g',alpha=0.25)
axes2.set_ylim(0, 1)
axes2.set_ylabel('Spatial Coverage (%)')
axes2.margins(x=0)
resfile=workpath+"/TS_nRMSE_nadirswot.png"
plt.savefig(resfile,bbox_inches="tight")    # save the figure
plt.close()             # close the figure


ymax_ = np.ceil(np.max([nrmse_Grad_OI_nadir, nrmse_Grad_OI_nadirswot, nrmse_Grad_Post_AnDA_nadir, nrmse_Grad_VE_DINEOF_nadir, nrmse_Grad_FP_ConvAE_nadir, nrmse_Grad_FP_GENN_nadir,\
                nrmse_Grad_Post_AnDA_nadirswot,nrmse_Grad_VE_DINEOF_nadirswot,nrmse_Grad_FP_ConvAE_nadirswot,nrmse_Grad_FP_GENN_nadirswot])*100)/100
# first axis with nRMSE time series
plt.plot(range(N),nrmse_Grad_OI_nadir,linestyle='solid',color='red',linewidth=2,label=r"$\nabla_{OI}$ (nadir)")
plt.plot(range(N),nrmse_Grad_Post_AnDA_nadir,linestyle='solid',color='seagreen',linewidth=1,label=r"$\nabla_{Post-AnDA}$ (nadir)")
plt.plot(range(N),nrmse_Grad_VE_DINEOF_nadir,linestyle='solid',color='steelblue',linewidth=1,markerSize=2,label=r"$\nabla_{VE-DINEOF}$ (nadir)")
plt.plot(range(N),nrmse_Grad_FP_ConvAE_nadir,linestyle='solid',color='mediumorchid',linewidth=1,markerSize=2,label=r"$\nabla_{FP-ConvAE}$ (nadir)")
plt.plot(range(N),nrmse_Grad_FP_GENN_nadir,linestyle='solid',color='darkorange',linewidth=1,markerSize=2,label=r"$\nabla_{FP-GENN}$ (nadir)")
plt.plot(range(N),nrmse_Grad_OI_nadirswot,linestyle='dashdot',color='red',linewidth=2,label=r"$\nabla_{OI}$ (nadir+swot)")
plt.plot(range(N),nrmse_Grad_Post_AnDA_nadirswot,linestyle='dashdot',color='seagreen',linewidth=1,label=r"$\nabla_{AnDA} (nadir+swot)$")
plt.plot(range(N),nrmse_Grad_VE_DINEOF_nadirswot,linestyle='dashdot',color='steelblue',linewidth=1,markerSize=2,label=r"$\nabla_{VE-DINEOF} (nadir+swot)$")
plt.plot(range(N),nrmse_Grad_FP_ConvAE_nadirswot,linestyle='dashdot',color='mediumorchid',linewidth=1,markerSize=2,label=r"$\nabla_{FP-ConvAE} (nadir+swot)$")
plt.plot(range(N),nrmse_Grad_FP_GENN_nadirswot,linestyle='dashdot',color='darkorange',linewidth=1,markerSize=2,label=r"$\nabla_{FP-GENN} (nadir+swot)$")
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
p2 = axes2.bar(range(N), nadswotcov-nadcov, width,bottom=nadcov,color='g',alpha=0.25)
axes2.set_ylim(0, 1)
axes2.set_ylabel('Spatial Coverage (%)')
axes2.margins(x=0)
resfile=workpath+"/TS_nRMSE_Grad.png"
plt.savefig(resfile,bbox_inches="tight")    # save the figure
plt.close()         	# close the figure

## Plot time series (nadir)
N = len(lday)
# first axis with nRMSE time series
plt.plot(range(N),nrmse_Grad_OI_nadir,linestyle='solid',color='red',linewidth=2,label=r"$\nabla_{OI}$ (nadir)")
plt.plot(range(N),nrmse_Grad_Post_AnDA_nadir,linestyle='solid',color='seagreen',linewidth=1,label=r"$\nabla_{Post-AnDA}$ (nadir)")
plt.plot(range(N),nrmse_Grad_VE_DINEOF_nadir,linestyle='solid',color='steelblue',linewidth=1,markerSize=2,label=r"$\nabla_{VE-DINEOF}$ (nadir)")
plt.plot(range(N),nrmse_Grad_FP_ConvAE_nadir,linestyle='solid',color='mediumorchid',linewidth=1,markerSize=2,label=r"$\nabla_{FP-ConvAE}$ (nadir)")
plt.plot(range(N),nrmse_Grad_FP_GENN_nadir,linestyle='solid',color='darkorange',linewidth=1,markerSize=2,label=r"$\nabla_{FP-GENN}$ (nadir)")
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
p2 = axes2.bar(range(N), nadswotcov-nadcov, width,bottom=nadcov,color='g',alpha=0.25)
axes2.set_ylim(0, 1)
axes2.set_ylabel('Spatial Coverage (%)')
axes2.margins(x=0)
resfile=workpath+"/TS_nRMSE_nadir_Grad.png"
plt.savefig(resfile,bbox_inches="tight")    # save the figure
plt.close()             # close the figure

## Plot time series (nadir/swot)
N = len(lday)
# first axis with nRMSE time series
plt.plot(range(N),nrmse_Grad_OI_nadirswot,linestyle='solid',color='red',linewidth=2,label=r"$\nabla_{OI}$ (nadir+swot)")
plt.plot(range(N),nrmse_Grad_Post_AnDA_nadirswot,linestyle='solid',color='seagreen',linewidth=1,label=r"$\nabla_{Post-AnDA}$ (nadir+swot)")
plt.plot(range(N),nrmse_Grad_VE_DINEOF_nadirswot,linestyle='solid',color='steelblue',linewidth=1,markerSize=2,label=r"$\nabla_{VE-DINEOF} (nadir+swot)$")
plt.plot(range(N),nrmse_Grad_FP_ConvAE_nadirswot,linestyle='solid',color='mediumorchid',linewidth=1,markerSize=2,label=r"$\nabla_{FP-ConvAE} (nadir+swot)$")
plt.plot(range(N),nrmse_Grad_FP_GENN_nadirswot,linestyle='solid',color='darkorange',linewidth=1,markerSize=2,label=r"$\nabla_{FP-GENN} (nadir+swot)$")
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
p2 = axes2.bar(range(N), nadswotcov-nadcov, width,bottom=nadcov,color='g',alpha=0.25)
axes2.set_ylim(0, 1)
axes2.set_ylabel('Spatial Coverage (%)')
axes2.margins(x=0)
resfile=workpath+"/TS_nRMSE_nadirswot_Grad.png"
plt.savefig(resfile,bbox_inches="tight")    # save the figure
plt.close()      



