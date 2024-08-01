import myfit

import glob
import os
import re
import sherpa.utils.err
from sherpa.astro.ui import *
from pychips import *

import pyntofale
import numpy as np

def myism(pars, x, y):
    return np.exp(-pars[0]*1e22*pyntofale.ionabs(x))

def load_data(maxorder):

    pha = 'output/_leg_m1.pha'
    bkg = 'output/_leg_m1_bkg.pha'

    set_default_id(1)
    load_pha(pha)#, use_errors=True)

    orders=range(1, maxorder+1)

    garfs=[ "output/m%d.arf" % (x) for x in orders ]
    rmfs=[ "output/m%d.rmf" % (x) for x in orders ]

    load_multi_arfs(garfs, orders)
    load_multi_rmfs(rmfs, orders)

    snr=30
    group_snr(snr, bkg_id=1)
    group_snr(snr)

    subtract()

    set_analysis('wave')

wlo = 12.39854 / 7.0
whi = 12.39854 / 0.2

load_data(10)
ignore(0., 10000.) ; notice(wlo, whi)

load_user_model(myism, 'ionabs')
add_user_pars('ionabs', ["nh"], [1.0])

#set_source(xstbabs.tbabs * xslogpar.lp)
#tbabs.nh = .0161
#freeze(tbabs.nh)

set_source(ionabs * xslogpar.lp)
ionabs.nh = .0161
freeze(ionabs.nh)

lp.alpha = 2
lp.beta = 0
lp.norm = 0.1

#fit()
lp.alpha = 2.3293
lp.beta = 0.0589843
lp.norm = 0.183262

title = 'Coadded Mkn 421: negative orders'

#myfit.plot_range(elow, ehigh, title)
#print_window( 'LEG_neg_full.pdf',  [ "clobber", True, 'fittopage', True, 'orientation', 'landscape'] )

#plot_data()

dplot = get_data_plot()
xdall = dplot.x
ydall = dplot.y
yerrdall = dplot.yerr

#plot_model()

mplot = get_model_plot()
xmall = 0.5 * (mplot.xlo + mplot.xhi)
ymall = mplot.y

#elow=0.25 ; ehigh=0.5 ; myfit.plot_range(elow, ehigh, title)
#print_window( 'LEG_neg.pdf',  [ "clobber", True, 'fittopage', True, 'orientation', 'landscape'] )

#elow=0.28 ; ehigh=0.31 ; myfit.plot_range(elow, ehigh, title)
#print_window( 'LEG_neg_carbon.pdf',  [ "clobber", True, 'fittopage', True, 'orientation', 'landscape'] )

load_data(1)
ignore(0., 10000.) ; notice(wlo, whi)

#elow=0.28 ; ehigh=0.31 ; myfit.plot_range(elow, ehigh, title)
#print_window( 'LEG_neg_carbon_1st_only.pdf',  [ "clobber", True, 'fittopage', True, 'orientation', 'landscape'] )

#plot_model()

mplot = get_model_plot()
xm1 = 0.5 * (mplot.xlo + mplot.xhi)
ym1 = mplot.y

# ugh, interpolation requires monotonically-increasing X values, and wavelengths are decreasing
xdall = xdall[::-1]
ydall = ydall[::-1]
yerrdall = yerrdall[::-1]

xmall = xmall[::-1]
ymall = ymall[::-1]

xm1 = xm1[::-1]
ym1 = ym1[::-1]

yd1 = ydall - np.interp(xdall, xmall, ymall-ym1)
yerrd1 = yerrdall

np.savetxt("foo.txt", np.column_stack((xdall, yd1, yerrd1, np.interp(xdall, xm1, ym1))))
