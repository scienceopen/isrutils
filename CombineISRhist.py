#!/usr/bin/env python3
"""
plots integrated ISR power over altitude range on top of HiST image stream

Michael Hirsch
"""
from __future__ import division, absolute_import
from six import integer_types
from pathlib2 import Path
from matplotlib.pyplot import show
from datetime import datetime
from pytz import UTC
#
from isrutils.common import boilerplateapi
from isrutils.plasmaline import readplasmaline#,plotplasmaline
from isrutils.summed import *
from GeoData.utilityfuncs import readNeoCMOS
#
heightkm=110.
epoch = datetime(1970,1,1,tzinfo=UTC)
#
def overlayisrhist(isrfn,odir,tlim,zlim,P):
    """
    1) read ISR plasma line
    2) read ISR long pulse
    3) sum over power from in NEIAL altitude range
    4) load HiST video
    5) register ISR to HST
    6) plot overlay joint data
    """
    assert isinstance(isrfn,Path)
    assert isinstance(zlim[0],(float,integer_types))

#%% (1) read ISR plasma line
#    plsum = sumplasmaline(isrfn,p.beamid,p.flim,tlim,zlim)
#    plotsumplasmaline(plsum)
    spec,freq = readplasmaline(isrfn,p.beamid,tlim)
    #plotplasmaline(spec,freq,isrfn)
#%% (2-3) read ISR long pulse
    lpsum,beamazel,isrlla = sumlongpulse(isrfn,p.beamid,tlim,zlim)
#%% (4) load optical data
    utlim = [(l-epoch).total_seconds() for l in tlim]

    hst = []; hstazel=[]; hstlla=[]; hstut=[]
    for f,a in zip(P.optfn,P.azelfn):
        opt, _, optazel, optlla, optut,_ = readNeoCMOS(P.optfn[0],P.azelfn[0],treq=utlim)
        hst.append(opt['optical']); hstazel.append(optazel)
        hstlla.append(optlla); hstut.append(optut)
#%% (5) transform magnetic zenith PFISR to HiST frame, assuming single altitude
    # now this happens inside do joint plot
#%% (6) plot joint
    dojointplot(lpsum,spec,freq,beamazel,hst,coordnames,optazel,optlla,isrlla,heightkm,
                utopt,utlim,P.makeplot)

if __name__ == '__main__':

    p,isrfn,odir,tlim = boilerplateapi()

    overlayisrhist(isrfn,odir,tlim,p.zlim,p)

    show()