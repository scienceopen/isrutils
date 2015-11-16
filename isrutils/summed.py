"""
summed measurements and plots
"""
from __future__ import division,absolute_import
from datetime import datetime
from pandas import Panel4D,DataFrame,Series
from numpy import absolute,nan,linspace,percentile,atleast_2d
from matplotlib.pyplot import figure,subplots,draw,pause
from matplotlib.cm import jet
from matplotlib.colors import LogNorm
#import matplotlib.animation as anim
#
from .plasmaline import readplasmaline
from .common import timeticks,findindex2Dsphere,timesync,projectisrhist
from .snrpower import readpower_samples
#from GeoData.plotting import plotazelscale

vidnorm = None #LogNorm()

#%% joint isr optical plot
def dojointplot(ds,spec,freq,beamazel,optical, optazel,optlla,isrlla,heightkm,utopt,utlim,makeplot):
    """
    f1,a1: radar   figure,axes
    f2,a2: optical figure,axes
    """
    assert isinstance(ds,(Series,DataFrame))
#%% setup master figure
    fig,axs = subplots(2,2,figsize=(8,12))
    axs = atleast_2d(axs)
#%% setup radar plot(s)
    a1 = axs[1,1]
    plotsumlongpulse(ds,a1)
    T = ds.index
    h1 = a1.axvline(nan,color='k',linestyle='--')
    t1 = a1.text(0.05,0.95,'time=',transform=a1.transAxes,va='top',ha='left')
#%% setup optical plot
    hi=[]; ht=[]; Iisr=[]; Iopt=[]
    for a,o,e,l in zip(axs[0,:],optical,optazel,optlla):
        optisrazel = projectisrhist(isrlla,beamazel,l,e,heightkm)

        clim = compclim(o,lower=10,upper=99.9)
        hi.append( a.imshow(o[0,...],origin='lower',interpolation='none',cmap='gray',
                       norm=vidnorm,vmin=clim[0],vmax=clim[1]) )
        a.set_axis_off()
        ht.append(a.set_title(''))

        azimg = e[:,1].reshape(o.shape[1:])
        elimg = e[:,2].reshape(o.shape[1:])
        #    plotazelscale(optical[0,...],azimg,elimg)

        #plot magnetic zenith beam
        br,bc = findindex2Dsphere(azimg,elimg,optisrazel['az'],optisrazel['el'])
        a.scatter(bc,br,s=500,marker='o',facecolors='none',edgecolor='red', alpha=0.5)

        a.autoscale(True,tight=True)


#%% setup plasma line plots
#    a2 = axs[1,0]; a3 = axs[1,1]
#    for a,c in zip((a2,a3),('down','up')):
#        plotplasmatime(spec,freq,t,fn,fig,a,tlim,vlim,c,makeplot,odir)

#%%
    #beam data, filled circle
#    s0 = a0.scatter(bc,br,s=500,alpha=0.5,linewidths=1.5,
#                    edgecolors=jet(linspace(ds.min(),ds.max())))

#%% time sync
        Is,Io = timesync(T,utopt,utlim)
        Iisr.append(Is); Iopt.append(Io)
#%% iterate
    for iisr,iopt0,iopt1 in zip(Iisr,Iopt[0],Iopt[1]):
#%% update isr plot
        t0isr=T[iisr]
        h1.set_xdata(t0isr)
        t1.set_text('isr: {}'.format(t0isr))
#%% update hist plot
        h0.set_data(optical[iopt,...])
        s0.set_array(ds[t0isr]) #magnetic zenith beam color (not changing color?)
        t0.set_text('optical: {}'.format(datetime.utcfromtimestamp(utopt[iopt])))
#%% anim
        if 'show' in makeplot:
            draw,pause(0.01)
        if 'png' in makeplot:
            fig.savefig('/tmp/joint_t{:05d}'.format(iopt),bbox_inches='tight',dpi=100)
#
#    def update(t):
#        h.set_xdata(t)
#        return h,
#
#    line_ani = anim.FuncAnimation(fig=f1, func=update, frames=T.size,
#                                   interval=50, blit=False)

def compclim(imgs,lower=0.5,upper=99.9,Nsamples=10):
    """
    inputs:
    images: Nframe x ypix x xpix grayscale image stack (have not tried with 4-D color)
    lower,upper: percentage (0,100)% to clip
    Nsamples: number of frames to test across the image stack (don't use too many for memory reasons)
    """
    sampind = linspace(0,imgs.shape[0],Nsamples,endpoint=False,dtype=int)

    return percentile(imgs[sampind,...],[lower,upper])

#%% dt3
def sumlongpulse(fn,beamid,tlim,zlim):
    snrsamp,azel,lla = readpower_samples(fn,beamid,tlim,zlim)
    assert isinstance(snrsamp,DataFrame)

    return snrsamp.sum(axis=0),azel,lla

def plotsumlongpulse(dsum,ax):
    assert isinstance(dsum,Series)
    if not ax:
        fg = figure()
        ax = fg.gca()

    dsum.plot(ax=ax)
    ax.set_ylabel('summed power')
    ax.set_xlabel('time [UTC]')
    ax.set_title('long pulse summed over altitude (200..350)km')

    ax.set_yscale('log')
    ax.grid(True)

    ax.xaxis.set_major_locator(timeticks(dsum.index[-1] - dsum.index[0]))
    return ax

#%% plasma line
def sumplasmaline(fn,beamid,flim,tlim,zlim):
    spec,freq = readplasmaline(fn,beamid,tlim)
    assert isinstance(spec,Panel4D)
    assert isinstance(flim[0],float)

    z = spec.major_axis
    specsum = DataFrame(index=spec.items,columns=spec.labels)

    zind = (zlim[0] <= z) & (z <= zlim[1])

    for s in spec:
        find = (flim[0] <= absolute(freq[s]/1e6)) & (absolute(freq[s]/1e6) < flim[1])
        specsum.loc[:,s] = spec.loc[:,:,zind,find].sum(axis=3).sum(axis=2)

    return specsum

def plotsumplasmaline(plsum):
    assert isinstance(plsum,DataFrame)

    fg = figure()
    ax = fg.gca()
    plsum.plot(ax=ax)
    ax.set_ylabel('summed power')
    ax.set_xlabel('time [UTC]')
    ax.set_title('plasma line summed over altitude (200..350)km and frequency (3.5..5.5)MHz')

    ax.xaxis.set_major_locator(timeticks(plsum.columns[-1]-plsum.columns[0]))
