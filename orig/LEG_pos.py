import misc

import glob
import os
import re
import sherpa.utils.err

base = 'data/mkn421/15483'

for d in ['data/mkn421/15483/fit/B0', 'data/mkn421/15483/fit/B0/binned']:
    try: os.mkdir(d)
    except: pass

pha2 = glob.glob(base+'/tg_reprocess/*_pha2.fits')
#pha2.extend( glob.glob(base+'/tg_reprocess/pha2_[0-9][0-9].fits'))

for pha2 in pha2:

    out='data/mkn421/15483/fit/B0'
    match = re.match( r'.*/pha2_(\d+).fits', pha2 )
    if match:
        out = out + '/binned/LEG_pos_' + match.groups()[0]
    else:
        out = out + '/LEG_pos'

    sdid=True

    if sdid:
        set_default_id(2)
        load_pha(pha2, use_errors = True)
    else:
        load_pha(2, pha2, use_errors = True)

    orders=range(1, 11)
    garfs=[ "data/mkn421/15483/tg_reprocess/hrcf15483N002LEG_%d_garf.fits" % (x) for x in orders ]
    rmfs=[ "data/mkn421/15483/tg_reprocess/LEG_%d.rmf" % (x) for x in orders ]

    load_multi_arfs(garfs, orders)
    load_multi_rmfs(rmfs, orders)

    snr=5
    group_snr(snr, bkg_id=1)
    group_snr(snr, bkg_id=2)
    group_snr(snr)

    if sdid:
        show_all( outfile='show_all_set_default_id', clobber=True )
    else:
        show_all( outfile='show_all_no_set_default_id', clobber=True )

    subtract()

    elow=0.5
    ehigh=8.0
    ignore(0., 10000.)
    notice(elow, ehigh)

    misc.source_def("mkn421", 0)


    fit()
    show_fit( outfile=out+'.params', clobber=True )

    misc.plot_range(elow, ehigh, "15483: TG\\_M = +1")

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

    if sdid:
        show_all( outfile='show_all_set_default_id', clobber=True )
    else:
        show_all( outfile='show_all_no_set_default_id', clobber=True )
