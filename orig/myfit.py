from sherpa.astro.ui import *
from pychips.extensions import *
from pychips.pychips_base import *
import os, errno

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise


def nH( src ):
    return {
        'mkn421'  : 0.0161,
        'pks2155' : 0.0142,
        }[src]

def source_def( src, band ):
    models = {
        0 : source_logpar,
    }
    models.get( band, source_powerlaw )(src)

def source_logpar( src ):
    set_source(xstbabs.tbabs * xslogpar.lp)
    tbabs.nh = nH( src )
    freeze(tbabs.nh)
    lp.alpha = 2
    lp.beta = 0
    lp.norm = 0.1

def source_powerlaw( src ):
    set_source(xstbabs.tbabs * xspowerlaw.pl)
    tbabs.nh = nH( src )
    freeze(tbabs.nh)
    pl.PhoIndex = 2.5
    pl.norm = 0.2

def plot_range(elow, ehigh, title, units='energy'):
    ignore()
    notice(elow, ehigh)
    plot(title, units)

def plot(title, units='kev'):

    set_ylog()

    plot_fit()
    split(2)
    adjust_grid_yrelsize(2,0.5)  # makes the second row 0.5 the height of the first
    current_plot('plot2')
    plot_ratio(overplot=True)

    bind_axes('plot1', 'ax1', 'plot2', 'ax1') # so that the X axes align

    if units == 'kev': set_plot_xlabel("Energy (keV)")
    else: set_plot_xlabel("Wavelength")
    set_plot_ylabel("Ratio (Data/Model)")

    current_plot('plot1')
    set_plot_title(title)
    hide_axis('ax1')
    if units == 'kev': log_scale(X_AXIS)

    current_plot('plot2')
    if units == 'kev': log_scale(X_AXIS)

    #plot_fit_delchi()
    #current_plot('plot2')
    #lin_scale(Y_AXIS)
