import myfit

import glob
import os
import re
import sherpa.utils.err
from sherpa.astro.ui import *
from pychips import *

base = 'data/mkn421/16430'

pha2 = glob.glob(base+'/tg_reprocess/*_pha2.fits')[0]

out='data/mkn421/16430/fit/B0'
match = re.match( r'.*/pha2_(\d+).fits', pha2 )
if match:
    out = out + '/binned/LEG_neg_' + match.groups()[0]
else:
    out = out + '/LEG_neg'

set_default_id(1)
load_pha(pha2, use_errors=True)

orders=range(1, 11)
garfs=[ "data/mkn421/16430/tg_reprocess/hrcf16430N002LEG_%d_garf.fits" % (-x) for x in orders ]
rmfs=[ "data/mkn421/16430/tg_reprocess/LEG_%d.rmf" % (-x) for x in orders ]

load_multi_arfs(garfs, orders)
load_multi_rmfs(rmfs, orders)

snr=20
group_snr(snr, bkg_id=1)
group_snr(snr, bkg_id=2)
group_snr(snr)

subtract()

wavlo=2.
wavhi=60.
set_analysis('wave')
ignore()
notice(wavlo, wavhi)

myfit.source_def("mkn421", 0)

fit()
show_fit( outfile=out+'.params', clobber=True )

flux = calc_energy_flux(wavlo, wavhi)
print("energy flux = %s\n" % flux)

# try:
#     conf()
#     show_conf()
# except sherpa.utils.err.EstErr:
#     pass

show_all()

myfit.plot_range(wavlo, wavhi, "16430: TG\_M = -1", units='wave')

print_window( '16430_gap_neg.pdf',  [ "clobber", True, 'fittopage', True, 'orientation', 'landscape'] )
