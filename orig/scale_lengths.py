import cPickle
import carmcmc as cm
import numpy as np
import sys
import argparse
import glob
import astropy.io.fits
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

def length_multiple(args):
    return { 'mkn421': 1e9, 'pks2155': 1e10 }[args.source]

def master(comm, data, sim_times, args):

    size = comm.Get_size()
    if size == 1:
        sys.stderr.write("there are no workers, silly\n")
        exit(1)

    status = MPI.Status()

    corr = np.arange(10*(args.maxscale-args.minscale)+1) * 10. / 100 + args.minscale
    scales = 1 + corr/100.

    return_means = np.zeros_like(sim_times)
    means = np.zeros((scales.size, sim_times.shape[0], sim_times.shape[1]))

    # get things rolling
    for i in xrange(min(size-1, scales.size)):
        comm.Send(scales[i], dest=i+1, tag=i)

    # wait for completed rows
    for i in xrange(size-1, scales.size):
        comm.Recv(return_means, source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
        means[status.Get_tag()] = return_means
        comm.Send(scales[i], dest=status.Get_source(), tag=i)

    # collect final results from each process
    for i in xrange(min(size-1, scales.size)):
        comm.Recv(return_means, source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
        means[status.Get_tag()] = return_means

    # send workers home
    for i in range(1, size):
        comm.Send(np.zeros(1)+np.nan, dest=i)

    # makes little difference which of these are used
    lengths = np.sqrt( (means[:,:,1:] - means[:,:,:-1] )**2).sum(axis=(1,2))
    lengths = np.sqrt( ( length_multiple(args) * (means[:,:,1:] - means[:,:,:-1]) )**2  + (sim_times[np.newaxis,:,1:] - sim_times[np.newaxis,:,:-1])**2).sum(axis=(1,2))

    for i in range(corr.size):
        pass
#        print scales[i], lengths[i]

    plot(data, sim_times, scales, means, args)

    return

def worker(comm, data, sim_times, args):
    status = MPI.Status()
    scales = np.zeros(1)

    while 1:

        comm.Recv(scales, source=0, tag=MPI.ANY_TAG, status=status)

        # all done if NaN was sent
        if scales != scales: return

        sim_means = get_sim_means(data, scales, sim_times, args)

        comm.Send(sim_means, dest=0, tag=status.Get_tag())

def get_sims(sample, times, scales, args):

    global npaths

    epoch = args.epoch
    source = args.source
    det1 = args.det1
    det2 = args.det2

    pickle_dir = 'data/pickle/B%d' % args.band
    pickle_file = '%s/%s_%s_%s_E%d_sims_%04.1f.bin' % (pickle_dir, source, det1, det2, epoch, 100*(scales-1))

    sims_shape = ( times.shape[0], npaths, times.shape[1] )
    try:
        f = open(pickle_file, 'r')
        sims = cPickle.load(f)
        f.close()

        if cmp(sims_shape, sims.shape):
            raise "The simulations must be recalculated"

    except:
        print(pickle_file)
        scale.mkdir_p(pickle_dir)
        sims = np.zeros(sims_shape)
        for i in range(sims.shape[0]):
            for j in range(sims.shape[1]):
                sims[i,j] = sample.simulate(times[i], bestfit='random')

        f = open(pickle_file, 'w')
        cPickle.dump(sims, f, 2)
        f.close()

    return sims

def get_samples(time, flux, err, scales, args):

    global nsamples

    epoch = args.epoch
    source = args.source
    det1 = args.det1
    det2 = args.det2

    pickle_dir = 'data/pickle/B%d' % args.band
    pickle_file = '%s/%s_%s_%s_E%d_samples_%04.1f.bin' % (pickle_dir, source, det1, det2, epoch, 100*(scales-1))

    try:
        f = open(pickle_file, 'r')
        sample = cPickle.load(f)
        f.close()

    except:
        print(pickle_file)
        scale.mkdir_p(pickle_dir)
        order_epochs, p, q = orders(args)
        if p.size-1 < epoch:
            raise "epoch requested greater than those for which orders have been chosen"
        model = cm.CarmaModel(time, np.log(flux), err/flux, p=int(p[epoch]), q=int(q[epoch]))
        sample = model.run_mcmc(nsamples)
        f = open(pickle_file, 'w')
        cPickle.dump(sample, f, 2)
        f.close()

    return sample

def get_sim_means(data, scales, sim_times, args):

    det1 = args.det1
    det2 = args.det2

    time = np.concatenate([data[det]['time'] for det in (det2, det1)])
    flux = np.concatenate((data[det2]['flux'], data[det1]['flux'] * scales))
    err = np.concatenate((data[det2]['err'], data[det1]['err'] * scales))
        
    samples = get_samples(time, flux, err, scales, args)

    sims = get_sims(samples, sim_times, scales, args)

    # a few scattered NaN are popping up in the simulations, and
    # they screw up everything afterwards
    sims = np.ma.masked_invalid(sims)

    # makes little difference whether we exponentiate before or after taking the mean
    mean = np.exp(sims.mean(axis=1))

    return mean

def orders(args):
    order_file='carma_orders_%s_%s_%s.rdb' % (args.source, args.det1, args.det2)
    order_epochs, p, q = np.loadtxt(order_file, skiprows=2, unpack=True )
    return order_epochs, p, q

def plot(data, sim_times, scales, means, args):

    band = args.band
    epoch = args.epoch
    source = args.source
    det1 = args.det1
    det2 = args.det2

    lengths = np.sqrt( ( length_multiple(args) * (means[:,:,1:] - means[:,:,:-1]) )**2  + (sim_times[np.newaxis,:,1:] - sim_times[np.newaxis,:,:-1])**2).sum(axis=(1,2))

    pdf = PdfPages('lengths_%s_%s_%s_B%d_E%d.pdf' % (source, det1, det2, band, epoch))

    if epoch == 0:
        title = 'all epochs'
    else:
        obsid = data[det1]['obsid'][0]
        file = glob.glob("data/%s/%d/tg_reprocess/*_evt2.fits" % (source, obsid))[0]
        hdulist = astropy.io.fits.open(file)
        title = hdulist[1].header['OBJECT'] + ': ' + hdulist[1].header['DATE-OBS'][0:10]

    fig = plt.figure(figsize = (11, 8.5))
    plt.plot(scales, lengths)
#    plt.xlim(1.0, 1.15)
    plt.title(title)
    plt.xlabel('scale factor')
    plt.ylabel('combined model path length: %s and scaled %s' % (det2, det1) )
    plt.tight_layout()
    pdf.savefig(fig)
    pdf.close()

    plt.close()

    if args.models:
        start = {}
        stop = {}

        for det in data:
            start[det], stop[det] = scale.get_gaps(data[det]['time'])
#            print det, start[det], stop[det]

        pdf = PdfPages('models_%s_%s_%s_B%d_E%d.pdf' % (source, det1, det2, band, epoch))
        fig = plt.figure(figsize = (11, 8.5))

        if epoch == 0:
            plot_dims = (2, int(math.ceil(start[det2].size/2.)))
        else:
            plot_dims = (1, 1)

        labels = { 'leg': 'ACIS-S/LETG', 
                   'meg': 'ACIS-S/MEG',
                   'heg': 'ACIS-S/HEG',
                   'hrc': 'HRC-S/LETG' }

        symbols = { 'leg': 'ko',
                    'meg': 'ro',
                    'heg': 'mo',
                    'hrc': 'bo' }

        for i in range(scales.size):

            for j in range(start[det2].size):

                row = int( j / plot_dims[1] )
                col = j % plot_dims[1]

                plt.subplot2grid(plot_dims, (row,col))

                for det in data:

                    starti = start[det][j]
                    stopi = stop[det][j] + 1

                    x = data[det]['time'][starti:stopi] - sim_times[j].min()
                    y = data[det]['flux'][starti:stopi]
                    yerr = data[det]['err'][starti:stopi]

                    if det == det1:
                        y = y * scales[i]
                        yerr = yerr * scales[i]

                    label = labels[det]
                    if det == det1:
                        label = 'scaled ' + label
                    plt.plot(x, y, symbols[det], label=label)
                    plt.errorbar(x, y, yerr, fmt=None)

                if row==plot_dims[0]-1: plt.xlabel(r'days')
                if col==0: plt.ylabel('%g - %g keV flux' % scale.band_range(band))

                plt.plot(sim_times[j]-sim_times[j].min(), means[i,j,:], 'k-')

                if j==0:
                    plt.title("scale = %.3f" % scales[i])
                    plt.legend(loc='lower left', fontsize=8)

            plt.tight_layout()
            pdf.savefig(fig)

        pdf.close()

        plt.close()

    return

def choose_order(data, args):
    det1 = args.det1
    det2 = args.det2

    scale = 1 + args.choose_order_scale/100.

    time = np.concatenate([data[det]['time'] for det in (det2, det1)])
    flux = np.concatenate((data[det2]['flux'], data[det1]['flux'] * scale))
    err = np.concatenate((data[det2]['err'], data[det1]['err'] * scale))

    model = cm.CarmaModel(time, np.log(flux), err/flux)
    model.choose_order(10, njobs=-1) # takes 30 minutes to run on legs

def main():

    global ntimes

    parser = argparse.ArgumentParser(
        description='Perform CARMA simulation models for two detectors\' combined fluxes, calculate length of model paths while scaling flux from the first detector.'
    )

    parser.add_argument('-e', '--epoch', type=int, default=0)
    parser.add_argument('-b', '--band', type=int, default=0)
    parser.add_argument('--minscale', type=float, default=2.)
    parser.add_argument('--maxscale', type=float, default=10.)
    parser.add_argument('-c', '--choose_order', action='store_true')
    parser.add_argument('--choose_order_scale', type=float, default=0.)
    parser.add_argument('--models', action='store_true')
    parser.add_argument('source')
    parser.add_argument('det1')
    parser.add_argument('det2')
    args = parser.parse_args()

    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()

    data = {}
    start = {}
    stop = {}
    det1 = args.det1
    det2 = args.det2

    for det in det1, det2:
        rdbfile = 'plots/%s_B%d_binned_%s.rdb' % (args.source, args.band, det)
        try:
            data[det] = scale.readflux( rdbfile )
        except IOError:
            if rank == 0: sys.stderr.write("could not read file '%s'" % rdbfile)
            exit(1)

    # FIXME - temporary fix for lone pks2155 acis data
    if det1 == 'hrc' or det2 == 'hrc':
        for det in data:
            mask = np.logical_or(data[det]['obsid'] < 8000, data[det]['obsid'] >= 9000)
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

#    exit(0)

    if args.choose_order:
        if rank: return
        choose_order(data, args)
        return

    sim_times = scale.get_sim_times(data, ntimes)

    if rank == 0:
        master(comm, data, sim_times, args)
    else:
        worker(comm, data, sim_times, args)

if __name__ == '__main__':
    main()
