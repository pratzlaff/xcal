import numpy as np
import argparse
import math
import astropy.io.fits
import glob

import scale

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from matplotlib import rc, rcParams
rc('text', usetex=True)

def main():

    parser = argparse.ArgumentParser(
        description='Calculate mean MEG/HEG ratio for the given obsids, using the spectral fluxes'
    )
    parser.add_argument('--combine', type=int, default=256)
    parser.add_argument('-p', '--pdf', help='Save plot to named file.')
    parser.add_argument('obsids', nargs='+')
    parser.add_argument('-t', '--title', help='Title for the plot.')
    args = parser.parse_args()

    spectra = { scale.HEG:{}, scale.MEG:{} }
    resp = { scale.HEG:{}, scale.MEG:{} }
    flux = { scale.HEG:{}, scale.MEG:{} }

    labels = {
        scale.HEG:{ 'pos':'HEG: TG\_M = +1', 'neg':'HEG: TG\_M = -1'},
        scale.MEG:{ 'pos':'MEG: TG\_M = +1', 'neg':'MEG: TG\_M = -1'}
        }

    w1, w2 = 1.8, 15.15
    for grating in spectra:
        combine = args.combine
        if grating == scale.HEG: combine *= 2
        spectra[grating] = scale.get_spectra(args.obsids, grating, combine)
        print spectra[grating]['pos']['lambda']
        for order in spectra[grating]:
            mask = np.logical_and(spectra[grating][order]['lambda']>=w1, spectra[grating][order]['lambda']<=w2)
            for s in 'lambda', 'flux', 'rate', 'flux_err':
                spectra[grating][order][s] = spectra[grating][order][s][mask]

    if args.pdf:
        pdf = PdfPages(args.pdf)

    f, axarr = plt.subplots(2, sharex=True, figsize=(8.5, 11))

    for grating in spectra:
        for order in spectra[grating]:
            x, y, yerr = spectra[grating][order]['lambda'], spectra[grating][order]['flux'], spectra[grating][order]['flux_err']
            axarr[0].plot(x, y, label=labels[grating][order])
            axarr[0].errorbar(x, y, yerr, ecolor='k', fmt=None)
    axarr[0].set_ylabel(r'$\textrm{Summed Flux (ergs/s/cm}^2\textrm{/\AA)}$')
    axarr[0].legend(loc='upper right', fontsize=8)
    if args.title:
        axarr[0].set_title(r'$\textrm{' + args.title + '}$')

    labels = { 'pos':'TG\_M = +1', 'neg':'TG\_M = -1'}
    ratio = {}
    ratio_err = {}
    for order in spectra[spectra.keys()[0]]:
        ratio[order] = spectra[scale.MEG][order]['flux'] / spectra[scale.HEG][order]['flux']
        ratio_err[order] = ratio[order] * np.sqrt((spectra[scale.MEG][order]['flux_err'] / spectra[scale.MEG][order]['flux'])**2 + (spectra[scale.HEG][order]['flux_err'] / spectra[scale.HEG][order]['flux'])**2)

        x, y, yerr = spectra[scale.MEG][order]['lambda'], ratio[order], ratio_err[order]
        axarr[1].plot(x, y, label=labels[order])
        axarr[1].errorbar(x, y, yerr, ecolor='k', fmt=None)
    axarr[1].legend(loc='upper left', fontsize=8)
    axarr[1].set_ylabel(r'$\textrm{MEG / HEG Ratio}$')
    axarr[1].set_xlabel(r'$\textrm{Wavelength (\AA)}$')

    f.tight_layout()
    f.subplots_adjust(hspace=0)

    if args.pdf:
        pdf.savefig(f)
        pdf.close()

    plt.show()

if __name__ == '__main__':
    main()
