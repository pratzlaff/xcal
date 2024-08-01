# see http://www.scholarpedia.org/article/Minimum_description_length
# and /data/fubar/SCAR/pro/stat/mdlpoly.pro

import cPickle
import carmcmc as cm
import numpy as np
from scipy.special import gammaln
import sys
import argparse
import glob
import pyfits
import math

import scale

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import matplotlib
#matplotlib.rc('text', usetex=True)

npaths=1000
nsamples=20000
ntimes=201

from mpi4py import MPI

# FIXME: only plots the first polynomial
def plot(data, sim_times, scales, orders, polys, args):

    band = args.band
    epoch = args.epoch
    source = args.source

    start = {}
    stop = {}

    for det in data:
        start[det], stop[det] = scale.get_gaps(data[det]['time'])

    pdf = PdfPages('mdl_%s_B%d_E%d_poly_%02d.pdf' % (source, band, epoch, orders[0]))
    fig = plt.figure(figsize = (11, 8.5))

    if epoch == 0:
        plot_dims = (2, int(math.ceil(start['leg'].size/2.)))
    else:
        plot_dims = (1, 1)

    labels = { 'leg': 'ACIS-S', 
               'hrc': 'scaled HRC-S' }
    symbols = { 'leg': 'ro',
                'hrc': 'bs' }

    for i in range(scales.size):

        for j in range(start['leg'].size):

            row = int( j / plot_dims[1] )
            col = j % plot_dims[1]

            plt.subplot2grid(plot_dims, (row,col))

            for det in data:

                starti = start[det][j]
                stopi = stop[det][j] + 1

                x = data[det]['time'][starti:stopi] - sim_times[j].min()
                y = data[det]['flux'][starti:stopi]
                yerr = data[det]['err'][starti:stopi]

                if det == 'hrc':
                    y = y * scales[i]
                    yerr = yerr * scales[i]

                plt.plot(x, y, symbols[det], label=labels[det])
                plt.errorbar(x, y, yerr, ecolor='k', fmt=None)

            if row==plot_dims[0]-1: plt.xlabel(r'days')
            if col==0: plt.ylabel('%g - %g keV flux' % scale.band_range(band))

            plt.plot(sim_times[j]-sim_times[j].min(), polys[0][i](sim_times[j]), 'k-')

            if j==0:
                plt.title("scale = %.3f" % scales[i])
                plt.legend(loc='lower left', fontsize=8)

        plt.tight_layout()
        pdf.savefig(fig)

    pdf.close()

    plt.close()

    return

def main():

    parser = argparse.ArgumentParser(description='Calculate the MDL for combined ACIS-S and scaled HRC-S LETG data, for a span of polynomial orders.')

    parser.add_argument('-e', '--epoch', type=int, default=0)
    parser.add_argument('-b', '--band', type=int, default=0)
    parser.add_argument('-o', '--order', type=int)
    parser.add_argument('--minscale', type=float, default=-5.)
    parser.add_argument('--maxscale', type=float, default=10.)
    parser.add_argument('--models', action='store_true')
    parser.add_argument('source')
    parser.add_argument('k1', type=int)
    parser.add_argument('k2', type=int)
    args = parser.parse_args()

    data = {}
    start = {}
    stop = {}

    for det in 'leg', 'hrc':
        rdbfile = 'plots/%s_B%d_binned_%s.rdb' % (args.source, args.band, det)
        try:
            data[det] = scale.readflux( rdbfile )
        except IOError:
            if rank == 0: sys.stderr.write("could not read file '%s'" % rdbfile)
            exit(1)

    # FIXME - temporary fix for lone pks2155 acis data
    for det in data:
        mask = data[det]['obsid'] != 8388
        for key in data[det]:
            data[det][key] = data[det][key][mask]

    # retain only the specified epoch
    epoch = args.epoch
    if epoch > 0:
        for det in data:
            start, stop = scale.get_gaps(data[det]['time'])
            if epoch > start.size:
                if rank == 0: sys.stderr.write("there are only %d epochs" % start.size)
                exit(1)
            for key in data[det]:
                data[det][key] = data[det][key][start[epoch-1]:stop[epoch-1]+1]

    # correction, in percent
    corr = np.arange(10*(args.maxscale-args.minscale)+1) * 10. / 100 + args.minscale
    scales = 1 + corr/100.

    mdl = np.zeros_like(scales)

    time = np.concatenate([data[det]['time'] for det in 'leg', 'hrc'])
    n = time.size

    orders = np.arange(args.k2-args.k1+1)+args.k1

    pdf = PdfPages('mdl_%s_B%d_E%d.pdf' % (args.source, args.band, args.epoch))
    fig = plt.figure(figsize = (11, 8.5))
    plot_dims=(2, 3)

    polys = []

    for j in xrange(orders.size):

        polys.append([])

        k = orders[j]

        for i in xrange(scales.size):
            flux = np.concatenate( (data['leg']['flux'], data['hrc']['flux'] * scales[i]) )
            err = np.concatenate(  (data['leg']['err'],  data['hrc']['err'] * scales[i] ) )

            coeff = np.polyfit(time, flux, k, w=1/err/err)
            p = np.poly1d(coeff)
            residuals = flux - p(time)

            polys[j].append(p)

            chisq = ((residuals/err)**2).sum()
            mdl[i] = chisq + k/2. * np.log(n) + gammaln(n+1) - gammaln(k+1) - gammaln(n-k+1)

        row = int( j / plot_dims[1] )
        col = j % plot_dims[1]

        min_mdl, min_scale = mdl.min(), scales[np.argmin(mdl)]

        plt.subplot2grid(plot_dims, (row,col))
        plt.plot(scales, mdl)
        plt.title('k=%d (%.4g, %.4g)' % (k,min_scale, min_mdl))
        plt.xlim(1.0, 1.1)

    plt.tight_layout()
    pdf.savefig(fig)
    pdf.close()

    if args.models:
        global ntimes
        sim_times = scale.get_sim_times(data, ntimes)
        plot(data, sim_times, scales, orders, polys, args)

    return

        
if __name__ == '__main__':
    main()
