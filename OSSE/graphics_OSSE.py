from tools_OSSE import *

## 1) Function for plotting nRMSE
def plot_nRMSE(list_data,labels_data,colors,symbols,lstyle,lwidth,lday,resfile,gradient=False):

    N   = len(lday)
    GT  = list_data[0]
    ## Compute spatial coverage  
    if any("Obs (nadir+swot)"==s for s in labels_data):
        spatial_coverage_nadirswot = []
        id_obs = np.where(labels_data=="Obs (nadir+swot)")[0][0]
        OBS = list_data[id_obs]
        for j in range(len(lday)):
            spatial_coverage_nadirswot.append(100*len(np.argwhere(np.isfinite(OBS[j].flatten())))/len(OBS[j].flatten()))
    if any("Obs (nadir)"==s for s in labels_data):
        spatial_coverage_nadir = []
        id_obs = np.where(labels_data=="Obs (nadir)")[0][0]
        OBS = list_data[id_obs]
        for j in range(len(lday)):
            spatial_coverage_nadir.append(100*len(np.argwhere(np.isfinite(OBS[j].flatten())))/len(OBS[j].flatten()))
        
    # Compute nRMSE time series
    nRMSE = []
    for i in range(len(labels_data[2:])):
        nRMSE_=[]
        meth_i=list_data[i+2]
        for j in range(len(lday)):
            if gradient == False:
                nRMSE_.append((np.sqrt(np.nanmean(((GT[j]-np.nanmean(GT[j]))-(meth_i[j]-np.nanmean(meth_i[j])))**2)))/np.nanstd(GT[j]))
            else:
                nRMSE_.append((np.sqrt(np.nanmean(((Gradient(GT[j],2)-np.nanmean(Gradient(GT[j],2)))-(Gradient(meth_i[j],2)-np.nanmean(Gradient(meth_i[j],2))))**2)))/np.nanstd(Gradient(GT[j],2)))
        nRMSE.append(nRMSE_)
  
    # plot nRMSE time series
    for i in range(len(labels_data[2:])):
        if gradient == False:
            plt.plot(range(N),nRMSE[i],linestyle=lstyle[i+2],color=colors[i+2],linewidth=lwidth[i+2],label=labels_data[i+2])
        else:
            plt.plot(range(N),nRMSE[i],linestyle=lstyle[i+2],color=colors[i+2],linewidth=lwidth[i+2],label=r"$\nabla_{"+str.split(labels_data[i+2])[0]+"}$ "+str.split(labels_data[i+2])[1])


    # add vertical bar to divide the 4 periods
    plt.axvline(x=19)
    plt.axvline(x=39)
    plt.axvline(x=59)
    # graphical options
    plt.ylim(0,np.ceil(np.max(nRMSE)*100)/100)
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
    if ( (any("Obs (nadir+swot)"==s for s in labels_data)) and any("Obs (nadir)"==s for s in labels_data) ):
        p1 = axes2.bar(range(N), spatial_coverage_nadir, width,color='r',alpha=0.25)
        p2 = axes2.bar(range(N), spatial_coverage_nadirswot-spatial_coverage_nadir, width,bottom=spatial_coverage_nadir,color='g',alpha=0.25)
    elif any("Obs (nadir+swot)"==s for s in labels_data):
        p1 = axes2.bar(range(N), spatial_coverage_nadirswot, width,color='g',alpha=0.25)
    else:
        p1 = axes2.bar(range(N), spatial_coverage_nadir, width,color='g',alpha=0.25)
    axes2.set_ylim(0, 100)
    axes2.set_ylabel('Spatial Coverage (%)')
    axes2.margins(x=0)
    plt.savefig(resfile,bbox_inches="tight")    # save the figure
    plt.close()                                 # close the figure

## 2) Function for plotting Signal-to-Noise ratio
def plot_SNR(list_data,labels_data,colors,symbols,lstyle,lwidth,lday,resssh,resfile):

    GT  = list_data[0]
    OBS = list_data[1]

    # Compute Signal-to-Noise ratio
    SNR = []
    for i in range(len(labels_data[2:])):
        f, Pf  = avg_err_raPsd2dv1(list_data[i],GT,resssh,True)
        wf     = 1./f
        SNR.append([wf, Pf])

    # plot Signal-to-Noise ratio
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for i in range(len(labels_data[2:])):
        wf, Pf = SNR[i]
        ax.plot(wf,Pf,linestyle=lstyle[i+2],color=colors[i+2],linewidth=lwidth[i+2],label=labels_data[i+2])
    ax.set_xlabel("Wavenumber", fontweight='bold')
    ax.set_ylabel("Signal-to-noise ratio", fontweight='bold')
    ax.set_xscale('log') ; ax.set_yscale('log')
    plt.legend(loc='best',prop=dict(size='small'),frameon=False)
    plt.xticks([50, 100, 200, 500, 1000], ["50km", "100km", "200km", "500km", "1000km"])
    ax.invert_xaxis()
    plt.grid(which='both', linestyle='--')
    plt.savefig(resfile) # save the figure
    plt.close()          # close the figure

## 3) Function for plotting Taylor diagrams
def Taylor_diagram(list_data,labels_data,colors,lstyle,resfile):

    # apply High-Pass Filter to visualize Taylor diagrams only for small scales
    HR = list_data[2]
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

    # create a dictionnary with data
    series={}
    for i in np.delete(np.arange(len(list_data)),1):
        series[labels_data[i]] = list_data[i].flatten()-lr
    Taylor_diag(series,np.delete(labels_data,1),\
            styles=np.delete(lstyle,1),\
            colors=np.delete(colors,1))
    plt.savefig(resfile)
    plt.close()


## 4) Display individual maps (SSH & Gradients)
def plot_maps(list_data,list_suffix,labels_data,list_day,date,extent,lon,lat,workpath):
    # vmin/vmax based on GT
    vmax = np.nanmax(np.abs(list_data[0]))
    vmin = -1.*vmax
    grad_vmax = np.nanmax(np.abs(Gradient(list_data[0],2)))
    grad_vmin = 0
    iday = np.where(list_day==date)[0][0]
    for i in range(len(list_data)):
        resfile = workpath+"/results_"+list_suffix[i]+'_'+date+".png"
        fig, ax = plt.subplots(1,1,figsize=(10,10),squeeze=False,
                          subplot_kw=dict(projection=ccrs.PlateCarree(central_longitude=0.0)))   
        plot(ax,0,0,lon,lat,list_data[i][iday],labels_data[i],\
             extent=extent,cmap="coolwarm",vmin=vmin,vmax=vmax)
        plt.savefig(resfile)       # save the figure
        plt.close()                # close the figure

        resfile = workpath+"/results_Grad_"+list_suffix[i]+'_'+date+".png"
        fig, ax = plt.subplots(1,1,figsize=(10,10),squeeze=False,
                          subplot_kw=dict(projection=ccrs.PlateCarree(central_longitude=0.0)))          
        if len(str.split(labels_data[i]))>1:
            plot(ax,0,0,lon,lat,Gradient(list_data[i][iday],2),r"$\nabla_{"+str.split(labels_data[i])[0]+"}$ "+str.split(labels_data[i])[1],\
                extent=extent,cmap="viridis",vmin=grad_vmin,vmax=grad_vmax)
        else:
            plot(ax,0,0,lon,lat,Gradient(list_data[i][iday],2),r"$\nabla_{"+labels_data[i]+"}$",\
                extent=extent,cmap="viridis",vmin=grad_vmin,vmax=grad_vmax)
        plt.savefig(resfile)       # save the figure
        plt.close()                # close the figure

## 5) Compute R/I/AE-scores
def RIAE_scores(list_data,labels_data,resfile,gradient=False):
    
    GT  = list_data[0]

    # apply High-Pass Filter to visualize Taylor diagrams only for small scales
    HR = list_data[2]
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

    def Iscore(mask1,gt,itrp):
        return 100*(1-np.nanmean(((mask1*gt-np.nanmean(mask1*gt))-(mask1*itrp-np.nanmean(mask1*itrp)))**2)/np.nanvar(mask1*gt))
    def Rscore(mask1,gt,itrp):
        return 100*(1-np.nanmean(((mask1*gt-np.nanmean(mask1*gt))-(mask1*itrp-np.nanmean(mask1*itrp)))**2)/np.nanvar(mask1*gt))
    def AEscore(gt,itrp):
        return 100*(1-np.nanmean(((gt-np.nanmean(gt))-(itrp-np.nanmean(itrp)))**2)/np.nanvar(gt))

    ## Compute mask nadir / swot  
    if any("Obs (nadir+swot)"==s for s in labels_data):
        spatial_coverage_nadir = []
        id_obs = np.where(labels_data=="Obs (nadir+swot)")[0][0]
        OBS = list_data[id_obs]
        mask1_nadirswot             = np.where(np.isnan(OBS.flatten()),np.nan,1)
        mask2_nadirswot             = np.where(np.isnan(OBS.flatten()),1,np.nan)
    if any("Obs (nadir)"==s for s in labels_data):
        spatial_coverage_nadir = []
        id_obs = np.where(labels_data=="Obs (nadir)")[0][0]
        OBS = list_data[id_obs]
        mask1_nadir                 = np.where(np.isnan(OBS).flatten(),np.nan,1)
        mask2_nadir                 = np.where(np.isnan(OBS.flatten()),1,np.nan)

    # Compute R/I/AE scores
    tab_scores = np.zeros((len(labels_data[2:]),3))
    for i in range(len(labels_data[2:])):
        meth_i=list_data[i+2]
        print(labels_data[i+2])
        if ("nadir+swot" in labels_data[i+2]):
            tab_scores[i,0] = Rscore(mask1_nadirswot,GT.flatten()-lr,meth_i.flatten()-lr)
            tab_scores[i,1] = Iscore(mask2_nadirswot,GT.flatten()-lr,meth_i.flatten()-lr)
        else:
            tab_scores[i,0] = Rscore(mask1_nadir,GT.flatten()-lr,meth_i.flatten()-lr)
            tab_scores[i,1] = Rscore(mask2_nadir,GT.flatten()-lr,meth_i.flatten()-lr)
        if ("rec" in labels_data[i+2]):
            tab_scores[i,2] = AEscore(GT.flatten()-lr,meth_i.flatten()-lr)
        else:
            tab_scores[i,2] = np.nan
    np.savetxt(fname=resfile,X=tab_scores,fmt='%2.2f')

