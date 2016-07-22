#!/usr/bin/env python
"""
plots integrated ISR power over altitude range on top of HiST image stream
"""
from isrutils.overlayISRopt import overlayisrhist

P={
'isrfn': '~/data/2013-04-14/isr/d0346834.dt3.h5',
'zlim': (200, 350),
'tlim': ['2013-04-14T08:54:10',
        '2013-04-14T08:54:20'],
'azelfn': '~/code/histfeas/precompute/hst0cal.h5',
'optfn':  '~/data/2013-04-14/hst/2013-04-14T8-54_hst0.h5',
'makeplot': [],
'beamid': 64157,
 'odir': 'out/2013-04-14'
}


overlayisrhist(P)
