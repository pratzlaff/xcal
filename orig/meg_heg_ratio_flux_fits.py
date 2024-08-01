import numpy as np
import argparse
import math
import astropy.io.fits
import glob

import scale

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from matplotlib import rc
#matplotlib.rc('text', usetex=True)

def plot_hist(ratios, title='', xlabel='', ylabel=''):
    plt.hist(ratios)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

def ratio_mean(x, xerr, y, yerr):
    ratios = x/y
    ratios_err = np.sqrt( ratios**2 * ( (xerr/x)**2 + (yerr/y)**2) )

    weights = 1/ratios_err/ratios_err

    weights /= np.sum(weights)

    mean = np.sum( ratios * weights )
    var = np.sum( weights**2 * ratios_err**2 )

    return mean, np.sqrt(var), ratios

def main():

    parser = argparse.ArgumentParser(
        description='Calculate mean MEG/HEG ratio for the given source, using the fitted fluxes'
    )

    parser.add_argument('source')
    args = parser.parse_args()

    title = { 'mkn421' : 'Mkn 421', 'pks2155' : 'PKS2155-304' }[args.source]
    pdf = PdfPages('%s_meg_heg_ratios.pdf' % args.source)
    fig = plt.figure(figsize = (11, 8.5))

    for band in (0, 3, 4):

        data = { }

        for det in 'meg', 'heg':
            rdbfile = 'plots/%s_B%d_binned_%s.rdb' % (args.source, band, det)
            try:
                data[det] = scale.readflux( rdbfile )
            except IOError:
                if rank == 0: sys.stderr.write("could not read file '%s'" % rdbfile)
                exit(1)

        start, stop = scale.get_gaps(data['meg']['time'])

        plot_dims = (2, int(math.ceil((start.size+1)/2.)))
        eband = scale.band_range(band)

        row, col = 0, 0
        mean, var, ratios = ratio_mean( data['meg']['flux'], data['meg']['err'], data['heg']['flux'], data['heg']['err'] )
        plt.subplot2grid(plot_dims, (row,col))
        plot_hist(ratios,
                  title=title + ': All epochs',
                  xlabel="%g - %g keV flux ratio: MEG / HEG" % (eband[0], eband[1]), ylabel='N')

        print(mean, var)

        for i in xrange(start.size):

            row = int( (i+1) / plot_dims[1] )
            col = (i+1) % plot_dims[1]

            tmp = { 'meg' : { }, 'heg' : { } }

            for det in tmp:
                for key in data[det]:
                    tmp[det][key] = data[det][key][start[i]:stop[i]+1]

            mean, var, ratios = ratio_mean( tmp['meg']['flux'], tmp['meg']['err'], tmp['heg']['flux'], tmp['heg']['err'] )

            print(mean, var)

            obsid = tmp['meg']['obsid'][0]
            file = glob.glob("data/%s/%d/tg_reprocess/*_evt2.fits" % (args.source, obsid))[0]
            hdulist = astropy.io.fits.open(file)
            date = hdulist[1].header['DATE-OBS'][0:10]

            plt.subplot2grid(plot_dims, (row,col))
            plot_hist(ratios, title=date)

        plt.tight_layout()
        pdf.savefig()

    pdf.close()
    plt.close()

if __name__ == '__main__':
    main()
