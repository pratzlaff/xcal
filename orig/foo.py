import sys, os, errno
sys.path=[os.environ['ASCDS_INSTALL'] + '/ots/lib/python' + ".".join(str(x) for x in sys.version_info[0:2]) + '/site-packages'] + sys.path 


import misc

import glob
from sherpa.astro.ui import *
from pychips import *

import re
import sherpa.utils.err

base = 'data/mkn421/10663'

for d in ['data/mkn421/10663/fit/B0', 'data/mkn421/10663/fit/B0/binned']:
    try: os.mkdir(d)
    except: pass

pha2 = glob.glob(base+'/tg_reprocess/*_pha2.fits')
pha2.extend( glob.glob(base+'/tg_reprocess/pha2_[0-9][0-9].fits')[0:2])

for pha2 in pha2:

    out='data/mkn421/10663/fit/B0'
    match = re.match( r'.*/pha2_(\d+).fits', pha2 )
    if match:
        out = out + '/binned/MEG_pos_' + match.groups()[0]
    else:
        out = out + '/MEG_pos'

    set_default_id(10)
    load_pha(pha2, use_errors=True)

    garf = "data/mkn421/10663/tg_reprocess/acisf10663N003MEG_1_garf.fits"
    rmf = "data/mkn421/10663/tg_reprocess/MEG_1.rmf"

    garf='../xcal/arfs/acis/10663_MEG_1_garf.fits'
    rmf='/data/legs/rpete/flight/rmfs/ACIS-S-MEG_1.rmf'
    olist=[1]

    load_multi_arfs((garf,), olist)
    load_multi_rmfs((rmf,), olist)

    # load_arf(garf)
    # load_arf(garf, bkg_id=1)
    # load_arf(garf, bkg_id=2)

    # load_rmf(rmf)

    snr=5
    group_snr(snr, bkg_id=1)
    group_snr(snr, bkg_id=2)
    group_snr(snr)

    subtract()

    elow=0.5
    ehigh=8.0
    ignore(0., 10000.)
    notice(elow, ehigh)

    misc.source_def("mkn421", 0)

    fit()
    #show_fit( outfile=out+'.params', clobber=True )

    flux = calc_energy_flux(elow, ehigh)
    #f=open( out+'.flux', 'w' )
    #f.write("energy flux = %s\n" % flux)
    #f.close()

    try:
	conf()
	#show_conf( outfile=out+'.conf', clobber=True )
    except sherpa.utils.err.EstErr:
	pass

    #show_all( outfile=out+'.show_all', clobber=True )

    #misc.plot_range(elow, ehigh, "10663: TG\_M = +1")

    #print_window( out+'.pdf',  [ "clobber", True, 'fittopage', True, 'orientation', 'landscape'] )
#    print_window( out+'.ps',  [ "clobber", True, 'fittopage', True, 'orientation', 'landscape'] )
#    print_window( out+'.png', [ "clobber", True ] )
