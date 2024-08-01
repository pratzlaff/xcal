import argparse

import myfit

import sherpa.utils.err
from sherpa.astro.ui import *
from pychips import *

import pyntofale
import numpy as np

def set_notice(wlo, whi):
    ignore()
    notice(wlo, whi)
    ignore(48, 58)

def mkplot(wlo, whi, filename=None):
    set_notice(wlo, whi)
    title = 'Coadded Mkn 421: negative orders'
    myfit.plot(title, units='ang')
    if filename is not None:
        print_window( filename,  [ "clobber", True, 'fittopage', True, 'orientation', 'landscape'] )

def myism(pars, x, y):

    nh = pars[0]
    cmult = pars[1]
    cii = pars[2]
    ciii = pars[3]
    oii = pars[4]
    oiii = pars[5]
    omult = pars[6]

    # Grevesse, N. \& Sauval, A.J.\ 1998, in Space Science Reviews,
    # 85, 161-174
    # CHIANTI says they are same as Grevesse, Noels, & Sauval 1996.
    abset = 'Grevesse & Sauval (1998)'
    abund = np.array([12.00, 10.93, 1.10, 1.40, 2.55, 8.52, 7.92, 8.83,
                      4.56, 8.08, 6.33, 7.58, 6.47, 7.55, 5.45, 7.33, 5.50, 6.40, 5.12,
                      6.36, 3.17, 5.02, 4.00, 5.67, 5.39, 7.50, 4.92, 6.25, 4.21, 4.60])
    abund = 10.**(abund-12.)
    abund[5] *= cmult
    abund[7] *= omult

    ionfrac = np.zeros((abund.size, abund.size+1))
    ionfrac[:,0] = 1

    ionfrac[5,0] = 1 - cii - ciii
    ionfrac[5,1] = cii
    ionfrac[5,2] = ciii

    ionfrac[7,0] = 1 - oii - oiii
    ionfrac[7,1] = oii
    ionfrac[7,2] = oiii

    return np.exp(-pars[0]*1e22*pyntofale.ionabs(x, abund=abund, ionfrac=ionfrac))

def load_data(maxorder, args):
    
    ddir = 'output/hrcs_qe' + args.qe

    pha = ddir + '/_leg_m1.pha'
    bkg = ddir + '/_leg_m1_bkg.pha'

    set_default_id(1)
    load_pha(pha)#, use_errors=True)

    orders=range(1, maxorder+1)

    garfs=[ ddir + '/m%d.arf' % (x) for x in orders ]
    rmfs=[ ddir + '/m%d.rmf' % (x) for x in orders ]

    load_multi_arfs(garfs, orders)
    load_multi_rmfs(rmfs, orders)

    snr=15
    group_snr(snr, bkg_id=1)
    group_snr(snr)

    subtract()

    set_analysis('wave')

def main():

    parser = argparse.ArgumentParser(
        description='Do some fitting of Mkn 421'
    )
    parser.add_argument('-i', '--ionabs',  action='store_true', default=False)
    parser.add_argument('qe')
    args = parser.parse_args()

    load_data(10, args)

    wlo=12.3985/7.0 ; whi=12.3985/0.2 ; set_notice(wlo, whi)

    if args.ionabs:
        load_user_model(myism, 'ionabs')
        add_user_pars('ionabs',
                      parnames = ['nh', 'cmult', 'cii', 'ciii', 'oii', 'oiii', 'omult'],
                      parvals = [1, 1, 0, 0, 0, 0, 1],
                      parmins = [0, 0.7, 0, 0, 0, 0, 0],
                      parmaxs = [10, 1.3, 0.5, 0.5, 0.5, 0.5, 1],
                      parfrozen = [False, True, True, True, True, True, True]
                      )
        set_source(ionabs * xslogpar.lp)
        nh = ionabs.nh

        ionabs.nh = .0161
        ionabs.nh.min = .01
        ionabs.nh.max = .025

    else:
        set_source(xstbabs.tbabs * xslogpar.lp)
        nh = tbabs.nh

        tbabs.nh = .0161
        tbabs.nh.min = .01
        tbabs.nh.max = .025
        # freeze(tbabs.nh)

    # nh = .0161
    # nh.min = .01
    # nh.max = .02

    lp.alpha = 2
    lp.beta = 0
    lp.norm = 0.1

    fit()

    postfix = ''
    if args.ionabs: postfix = '_ionabs'

    plotdir = './plots/hrcs_qe' + args.qe
    myfit.mkdir_p(plotdir)

    mkplot(5, 65, plotdir+'/LEG_neg_full' + postfix + '.pdf')
    mkplot(25, 47, plotdir+'/LEG_neg' + postfix + '.pdf')
    mkplot(41, 44, plotdir+'/LEG_neg_C' + postfix + '.pdf')
    mkplot(22, 24, plotdir+'/LEG_neg_O' + postfix + '.pdf')

    load_data(1, args)
    mkplot(41, 44, plotdir+'/LEG_neg_C_1st_only' + postfix + '.pdf')

    if args.ionabs:

        load_data(10, args)

        set_notice(41, 44)
        thaw(ionabs.cmult, ionabs.cii, ionabs.ciii)
        freeze(lp.alpha, lp.beta, lp.norm, ionabs.nh)
        fit()
        mkplot(41, 44, plotdir+'/LEG_neg_C' + postfix + '_fit.pdf')

        set_notice(22, 24)
        freeze(ionabs.cmult, ionabs.cii, ionabs.ciii)
        thaw(ionabs.oii, ionabs.oiii)
        fit()
        mkplot(22, 24, plotdir+'/LEG_neg_O' + postfix + '_fit.pdf')

    pass

if __name__ == '__main__':
    main()
