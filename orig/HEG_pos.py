import misc

import glob
import os
import re
import sherpa.utils.err
import numpy

base = 'data/mkn421/13105'

for d in ['data/mkn421/13105/fit/B3', 'data/mkn421/13105/fit/B3/binned']:
    try: os.mkdir(d)
    except: pass

#pha2 = glob.glob(base+'/tg_reprocess/*_pha2.fits')
pha2 = glob.glob(base+'/tg_reprocess/pha2_00.fits')

for pha2 in pha2:

    out='data/mkn421/13105/fit/B3'
    match = re.match( r'.*/pha2_(\d+).fits', pha2 )
    if match:
        out = out + '/binned/HEG_pos_' + match.groups()[0]
    else:
        out = out + '/HEG_pos'

    set_default_id(4)
    load_pha(pha2, use_errors = True)

    load_arf("data/mkn421/13105/tg_reprocess/acisf13105N003HEG_1_garf.fits");
    load_arf("data/mkn421/13105/tg_reprocess/acisf13105N003HEG_1_garf.fits", bkg_id=1);
    load_arf("data/mkn421/13105/tg_reprocess/acisf13105N003HEG_1_garf.fits", bkg_id=2);

    load_rmf("data/mkn421/13105/tg_reprocess/HEG_1.rmf");

    snr=5
    group_snr(snr, bkg_id=1)
    group_snr(snr, bkg_id=2)
    group_snr(snr)

    subtract()

    elow=0.85
    ehigh=1.5
    ignore(0., 10000.)
    notice(elow, ehigh)

    misc.source_def("mkn421", 3)

    fit()

    # a = sample_energy_flux(elow, ehigh, num=1000)
    # numpy.ndarray.sort( a[:,0] )
    # mean = a[:,0].mean()
    # std = a[:,0].std()
    # print("energy flux = %g (%g)" % (mean, std))

    # print("lopped distribution energy flux is [%g, %g]" % (a[160:161,0], a[840:841,0]))

    conf()


#     show_fit( outfile=out+'.params', clobber=True )

#     flux = calc_energy_flux(elow, ehigh)
#     f=open( out+'.flux', 'w' )
#     f.write("energy flux = %s\n" % flux)
#     f.close()

#     try:
# 	conf()
# 	show_conf( outfile=out+'.conf', clobber=True )
#     except sherpa.utils.err.EstErr:
# 	pass
#     misc.plot_range(elow, ehigh, "13105: TG\_M = +1")

#     print_window( out+'.pdf',  [ "clobber", True, 'fittopage', True, 'orientation', 'landscape'] )
# #    print_window( out+'.ps',  [ "clobber", True, 'fittopage', True, 'orientation', 'landscape'] )
# #    print_window( out+'.png', [ "clobber", True ] )

