import misc

import glob
import os
import re
import sherpa.utils.err

base = 'data/pks2155/1705'

for d in ['data/pks2155/1705/fit/B0', 'data/pks2155/1705/fit/B0/binned']:
    try: os.mkdir(d)
    except: pass

pha2 = glob.glob(base+'/tg_reprocess/*_pha2.fits')
#pha2.extend( glob.glob(base+'/tg_reprocess/pha2_[0-9][0-9].fits'))

for pha2 in pha2:

    out='data/pks2155/1705/fit/B0'
    match = re.match( r'.*/pha2_(\d+).fits', pha2 )
    if match:
        out = out + '/binned/MEG_pos_' + match.groups()[0]
    else:
        out = out + '/MEG_pos'

    set_default_id(10)
    load_pha(pha2, use_errors = True)

    load_arf("data/pks2155/1705/tg_reprocess/acisf01705N005MEG_1_garf.fits");
    load_arf("data/pks2155/1705/tg_reprocess/acisf01705N005MEG_1_garf.fits", bkg_id=1);
    load_arf("data/pks2155/1705/tg_reprocess/acisf01705N005MEG_1_garf.fits", bkg_id=2);

    load_rmf("data/pks2155/1705/tg_reprocess/MEG_1.rmf");

    snr=5
    group_snr(snr, bkg_id=1)
    group_snr(snr, bkg_id=2)
    group_snr(snr)

    subtract()

    elow=0.5
    ehigh=8.0
    ignore(0., 10000.)
    notice(elow, ehigh)

    misc.source_def("pks2155", 0)

    fit()
    show_fit( outfile=out+'.params', clobber=True )

    misc.plot_range(elow, ehigh, "1705: TG\_M = +1")

    print_window( out+'.pdf',  [ "clobber", True, 'fittopage', True, 'orientation', 'landscape'] )
#    print_window( out+'.ps',  [ "clobber", True, 'fittopage', True, 'orientation', 'landscape'] )
#    print_window( out+'.png', [ "clobber", True ] )

    flux = calc_energy_flux(elow, ehigh)
    f=open( out+'.flux', 'w' )
    f.write("energy flux = %s\n" % flux)
    f.close()

    try:
	conf()
	show_conf( outfile=out+'.conf', clobber=True )
    except sherpa.utils.err.EstErr:
	pass
