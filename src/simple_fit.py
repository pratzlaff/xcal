import glob

from sherpa.astro.ui import *

# from https://cxc.cfa.harvard.edu/sherpa/threads/grating_hrcsletg/index.html#response

import numpy as np
from sherpa.models.model import CompositeModel
from sherpa.astro.instrument import MultiResponseSumModel

def startup_monkey(self, cache):
    pha = self.pha
    if np.iterable(pha.mask):
        pha.notice_response(True)
    self.channel = pha.get_noticed_channels()
    self.mask = pha.get_mask()
    self._get_noticed_energy_list()
    CompositeModel.startup(self, cache)

MultiResponseSumModel.startup = startup_monkey

obsid=10665
src='mkn421'

row=2
pha2=glob.glob('/data/legs/rpete/flight/xcal/data/{}/{}/tg_reprocess/*_pha2.fits'.format(src, obsid))[0]

set_default_id(row)
load_pha(pha2, use_errors=True)

group_counts(1)

set_stat('wstat')

if False:
    garf='/data/legs/rpete/flight/xcal/arfs/qe_N0014_qeu_N0013/{}_LEG_{}_garf.fits'.format(obsid, 1)
    rmf='/data/legs/rpete/flight/rmfs/HRC-S-LEG_{}.rmf'.format(1)
    load_arf(garf)
    load_rmf(rmf)
else:
    olist=range(1,11)
    garfs=['/data/legs/rpete/flight/xcal/arfs/qe_N0014_qeu_N0013/{}_LEG_{}_garf.fits'.format(obsid, o) for o in olist]
    rmfs=['/data/legs/rpete/flight/rmfs/HRC-S-LEG_{}.rmf'.format(o) for o in olist]
    load_multi_arfs(garfs, olist)
    load_multi_rmfs(rmfs, olist)

elo, ehi = 0.54, 0.85 # kev
wlo, whi = 12.39854/ehi, 12.39854/elo # AA

set_analysis('wave')
ignore()
notice(wlo, whi)

set_xsabund("wilm")
set_source(xstbabs.tbabs * xspowerlaw.pl)
tbabs.nh = 0.0161
freeze(tbabs.nh)
set_par(pl.PhoIndex, val=2.5, min=-2, max=9)
set_par(pl.norm, val=0.2, min=0, max=1e24)

fit()
