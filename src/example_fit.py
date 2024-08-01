# Sherpa wstat throws a lot of background warnings
import warnings
warnings.simplefilter('ignore', UserWarning)

title='10665: LEG, TG\_M=+1'
olist=[1,2,3,4,5,6,7,8,9,10]
row=2
pha2='/data/legs/rpete/flight/xcal/data/mkn421/10665/tg_reprocess/hrcf10665N003_pha2.fits'
garfs=['/data/legs/rpete/flight/xcal/arfs/qe_N0014_qeu_N0011/10665_LEG_{}_garf.fits'.format(o) for o in olist]
rmfs=['/data/legs/rpete/flight/rmfs/HRC-S-LEG_{}.rmf'.format(o) for o in olist]

title='9712: MEG, TG\_M=+1'
olist=[1]
row=10
pha2='/data/legs/rpete/flight/xcal/data/pks2155/9712/tg_reprocess/acisf09712N003_pha2.fits'
garfs=['/data/legs/rpete/flight/xcal/arfs/acis/9712_MEG_{}_garf.fits'.format(o) for o in olist]
rmfs=['/data/legs/rpete/flight/rmfs/ACIS-S-MEG_{}.rmf'.format(o) for o in olist]

set_default_id(row)
load_pha(pha2, use_errors=True)

group_counts(1)

set_stat('wstat')

load_multi_arfs(garfs, olist)
load_multi_rmfs(rmfs, olist)

elo, ehi = 0.54, 0.85 # kev
wlo, whi = 12.39854/ehi, 12.39854/elo # AA

set_analysis('wave')
ignore()
notice(wlo, whi)

set_xsabund("wilm")
set_source(xstbabs.tbabs * xspowerlaw.pl)
tbabs.nh = 0.0142
freeze(tbabs.nh)
set_par(pl.PhoIndex, val=2.5, min=-2, max=9)
set_par(pl.norm, val=0.2, min=0, max=1e24)

fit()
#show_fit()

flux = calc_energy_flux(wlo, whi)
print('calc_energy_flux({:.3g}, {:.3g}) = {:g}'.format(wlo, whi, flux))

snr=5
group_snr(snr)
show_all()

conf()

if False:
    s = sample_flux(pl, wlo, whi, num=1000)

    fluxes = np.sort(s[2][:,0])
    np.save('fluxes.npy', fluxes)
    fluxes = np.load('fluxes.npy')

    offset=int(fluxes.size*(1-.6827)/2)
    s[0]
    print('flues[offset], fluxes[-offset] = {}, {}'.format(fluxes[offset], fluxes[-offset]))
    print('flues[offset+1], fluxes[-offset-4] = {}, {}'.format(fluxes[offset+1], fluxes[-offset-4]))

def dual_plot(title):
    units='wavelength' # or 'energy'

    set_ylog()

    plot_fit()
    split(2)
    adjust_grid_yrelsize(2,0.5)  # makes the second row 0.5 the height of the first
    current_plot('plot2')
    plot_ratio(overplot=True)

    bind_axes('plot1', 'ax1', 'plot2', 'ax1') # so that the X axes align

    if units == 'energy': set_plot_xlabel("Energy (keV)")
    else: set_plot_xlabel("Wavelength")
    set_plot_ylabel("Ratio (Data/Model)")

    current_plot('plot1')
    set_plot_title(title)
    hide_axis('ax1')
    if units == 'energy': log_scale(X_AXIS)

    current_plot('plot2')
    if units == 'energy': log_scale(X_AXIS)

    #plot_fit_delchi()
    #current_plot('plot2')
    #lin_scale(Y_AXIS)

dual_plot(title)
ungroup()
group_counts(1)

