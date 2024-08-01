import argparse
import astropy.io.fits
import cPickle
import glob
import math
from mpi4py import MPI
import numpy as np
import sys

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from matplotlib import rc, rcParams
#matplotlib.rc('text', usetex=True)
#rcParams.update({'font.size': 14})

import carmcmc as cm

import scale as scale_module
import xcal

rank=None

npaths=1000
nsamples=20000
ntimes=201

def length_multiple(source):
    return { 'mkn421': 1e9, 'pks2155': 1e10 }[source]

def master(comm, data, sim_days, args, hdr):

    size = comm.Get_size()
    if size == 1:
        sys.stderr.write("there are no workers, silly\n")
        exit(1)

    status = MPI.Status()

    scales = 1/100. * (np.arange((args.maxscale-args.minscale)/args.increment+1) * args.increment + args.minscale)
    #scales = 1 + corr/100.

    return_vals = np.zeros((sim_days.shape[0], sim_days.shape[1]+1))
    means = np.zeros((scales.size, sim_days.shape[0], sim_days.shape[1]))
    chi2 = np.zeros((scales.size, sim_days.shape[0]))

    # get things rolling
    for i in xrange(min(size-1, scales.size)):
        comm.Send(scales[i], dest=i+1, tag=i)

    # wait for completed rows
    for i in xrange(size-1, scales.size):
        comm.Recv(return_vals, source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
        means[status.Get_tag()] = return_vals[:,:-1]
        chi2[status.Get_tag()] = return_vals[:,-1]
        comm.Send(scales[i], dest=status.Get_source(), tag=i)

    # collect final results from each process
    for i in xrange(min(size-1, scales.size)):
        comm.Recv(return_vals, source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
        means[status.Get_tag()] = return_vals[:,:-1]
        chi2[status.Get_tag()] = return_vals[:,-1]

    # send workers home
    for i in range(1, size):
        comm.Send(np.zeros(1)+np.nan, dest=i)

#     # makes little difference which of these are used
#     lengths = np.sqrt( (means[:,:,1:] - means[:,:,:-1] )**2).sum(axis=(1,2))
#     lengths = np.sqrt( ( length_multiple(hdr['source']) * (means[:,:,1:] - means[:,:,:-1]) )**2  + (sim_days[np.newaxis,:,1:] - sim_days[np.newaxis,:,:-1])**2).sum(axis=(1,2))

#     for i in range(scales.size):
#         pass
# #        print scales[i], lengths[i]

    plot(data, sim_days, scales, means, chi2, args, hdr)

    return

def worker(comm, data, sim_days, args, hdr):
    status = MPI.Status()
    scales = np.zeros(1)

    while 1:

        comm.Recv(scales, source=0, tag=MPI.ANY_TAG, status=status)

        # all done if NaN was sent
        if scales != scales: return

        sim_means, sim_chi2 = get_sim_means(data, scales[0], sim_days, args, hdr)
        return_vals = np.zeros(sim_means.size+1)
        return_vals[:-1] = sim_means
        return_vals[-1] = sim_chi2

        comm.Send(return_vals, dest=0, tag=status.Get_tag())

def get_sims(sample, days, scale, args, hdr):

    global npaths

    epoch = args.epoch
    source = hdr['source']
    i1 = args.inst1
    i2 = args.inst2

    pickle_file = '{}/{}_{}{}_B{}_E{:02d}_sims_{:05.1f}.bin'.format(args.pdir, hdr['source'], args.inst1, args.inst2, hdr['band'], epoch, 100*scale)

    sims_shape = ( days.shape[0], npaths, days.shape[1] )
    try:
        f = open(pickle_file, 'rb')
        sims = cPickle.load(f)
        f.close()

        if cmp(sims_shape, sims.shape):
            raise "The simulations must be recalculated"

    except:
        print(pickle_file)
        xcal.mkdir_p(args.pdir)
        sims = np.zeros(sims_shape)
        for i in range(sims.shape[0]):
            for j in range(sims.shape[1]):
                sims[i,j] = sample.simulate(days[i], bestfit='random')

        f = open(pickle_file, 'wb')
        cPickle.dump(sims, f, 2)
        f.close()

    # a few scattered NaN are popping up in the simulations, and
    # they screw up everything afterwards
    return np.ma.masked_invalid(sims)

def get_samples(days, flux, err, scale, args, hdr):
    global nsamples

    pickle_file = '{}/{}_{}{}_B{}_E{:02d}_samples_{:05.1f}.bin'.format(args.pdir, hdr['source'], args.inst1, args.inst2, hdr['band'], args.epoch, 100*scale)

    try:
        f = open(pickle_file, 'rb')
        sample = cPickle.load(f)
        f.close()

    except:
        print(pickle_file)
        xcal.mkdir_p(args.pdir)
        bands, epochs, p, q = orders(args, hdr)
        i = np.where((epochs==args.epoch) & (bands==int(hdr['band'])))[0][0]
        if p.size-1 < args.epoch:
            raise "epoch requested greater than those for which orders have been chosen"
        p, q = p[i], q[i]
        # seeing weird results with p,q=1,0
        if p==1 and q==0:
            p=2

        if args.nolog:
            model = cm.CarmaModel(days, flux, err, p=int(p), q=int(q))
        else:
            model = cm.CarmaModel(days, np.log(flux), err/flux, p=int(p), q=int(q))
        sample = model.run_mcmc(nsamples)
        f = open(pickle_file, 'wb')
        cPickle.dump(sample, f, 2)
        f.close()

    return sample

def get_sim_means(data, scale, sim_days, args, hdr):

    i1 = args.inst1
    i2 = args.inst2

    days = np.append(data[i2]['days'], data[i1]['days'])
    flux = np.append(data[i2]['flux'], data[i1]['flux']*scale)
    err = np.append(data[i2]['err'], data[i1]['err']*scale)
        
    samples = get_samples(days, flux, err, scale, args, hdr)

    sims = get_sims(samples, sim_days, scale, args, hdr)

    mean = sims.mean(axis=1)
    # makes little difference whether we exponentiate before or after taking the mean
    if not args.nolog:
        mean = np.exp(mean)

    chi2 = np.zeros(mean.shape[0])
    for i in range(chi2.size):
        mean_interp = np.interp(days, sim_days[i], mean[i])
        chi2[i] = np.sum( ((mean_interp - flux) / err)**2)

    return mean, chi2

def orders(args, hdr):
    order_file='./carma_orders/{}_{}{}.rdb'.format(hdr['source'], args.inst1, args.inst2)
    bands, epochs, p, q = np.loadtxt(order_file, skiprows=2, unpack=True )

    return bands, epochs, p, q

def plot(data, sim_days, scales, means, chi2, args, hdr):
    band = int(hdr['band'])
    epoch = args.epoch
    source = hdr['source']
    i1 = args.inst1
    i2 = args.inst2

    lengths = np.sqrt( ( length_multiple(hdr['source']) * (means[:,:,1:] - means[:,:,:-1]) )**2  + (sim_days[np.newaxis,:,1:] - sim_days[np.newaxis,:,:-1])**2).sum(axis=(1,2))

    xcal.mkdir_p(args.pdfdir)
    pdf = PdfPages('{}/lengths_{}_{}{}_B{}_E{:02d}.pdf'.format(args.pdfdir, source, i1, i2, band, epoch))

    if epoch == 0:
        title = 'all epochs'
    else:
        obsid = data[i1]['obsid'][0]
        file = glob.glob("data/%s/%d/tg_reprocess/*_evt2.fits" % (source, obsid))[0]
        hdulist = astropy.io.fits.open(file)
        title = hdulist[1].header['OBJECT'] + ': ' + hdulist[1].header['DATE-OBS'][0:10] + ', {} and Scaled {}'.format(i2, i1)


    if args.chi2:
        fig = plt.figure(figsize = (8.5, 11))
        plt.subplot(211)
    else:
        fig = plt.figure(figsize = (11, 8.5))
        plt.xlabel('Scale Factor')

    plt.plot(scales, lengths)
    plt.title(title)
    plt.ylabel('Model Path Length')

    if args.chi2:
        plt.subplot(212)
        plt.xlabel('Scale Factor')
        plt.plot(scales, chi2)
        plt.xlabel('Scale Factor')
        plt.ylabel('$\chi^2$')

    plt.tight_layout()
    pdf.savefig(fig)
    pdf.close()

    plt.close()

    if args.models:
        groups = { inst : xcal.group_days(data[inst]['days']) for inst in data }

        pdf = PdfPages('{}/models_{}_{}{}_B{}_E{:02d}.pdf'.format(args.pdfdir, source, i1, i2, band, epoch))
        fig = plt.figure(figsize = (11, 8.5))

        ngroups = len(groups[list(groups.keys())[0]])

        if epoch == 0:
            plot_dims = (2, int(math.ceil(ngroups/2.)))
        else:
            plot_dims = (1, 1)

        labels = { 'l': 'ACIS/LEG', 
                   'm': 'ACIS/MEG',
                   'h': 'ACIS/HEG',
                   's': 'HRC/LEG' }

        symbols = { 'l': 'ko',
                    'm': 'ro',
                    'h': 'mo',
                    's': 'bo' }

        for i in range(scales.size):

            for j in range(ngroups):
                
                row = int( j / plot_dims[1] )
                col = j % plot_dims[1]

                plt.subplot2grid(plot_dims, (row,col))

                for inst in data:
                    g = groups[inst][j]

                    x = data[inst]['days'][g] - sim_days[j].min()
                    y = data[inst]['flux'][g]
                    yerr = data[inst]['err'][g]

                    label = labels[inst]

                    if inst == i1:
                        y = y * scales[i]
                        yerr = yerr * scales[i]
                        label = 'Scaled ' + label

                    #plt.plot(x, y, symbols[inst], label=label)
                    plt.errorbar(x, y, yerr, fmt=symbols[inst], label=label)

                if row==plot_dims[0]-1:
                    plt.xlabel(r'Days')
                if col==0:
                    elo, ehi = xcal.band_energy_limits(band)
                    plt.ylabel('{:.2g} - {:.2g} keV Flux'.format(elo, ehi))

                plt.plot(sim_days[j]-sim_days[j].min(), means[i,j,:], 'k-')

                if j==0:
                    plt.title("Scale = %.3f" % scales[i])
                    plt.legend(loc='lower left', fontsize=8)

            plt.tight_layout()
            pdf.savefig(fig)

        pdf.close()

        plt.close()

    return

def choose_order(data, args):
    i1 = args.inst1
    i2 = args.inst2

    scale = args.choose_order_scale/100.

    days = np.append(data[i2]['days'], data[i1]['days'])
    flux = np.append(data[i2]['flux'], data[i1]['flux'] * scale)
    err = np.append(data[i2]['err'], data[i1]['err'] * scale)

    if args.nolog:
        model = cm.CarmaModel(days, flux, err)
    else:
        model = cm.CarmaModel(days, np.log(flux), err/flux)
    model.choose_order(5, njobs=-1) # takes 30 minutes to run on legs

def get_sim_days(data, n):

    days = np.concatenate(([data[inst]['days'] for inst in data]))
    days = np.sort(days)

    groups = xcal.group_days(days)
    sim_days = np.zeros((len(groups), n))

    for i, g in enumerate(groups):
        span = days[g[-1]]-days[g[0]]
        sim_days[i]= np.arange(n) * span / n + days[g[0]]

    return sim_days

def main():

    global rank
    global ntimes

    parser = argparse.ArgumentParser(
        description='Perform CARMA simulation models for two detectors\' combined fluxes, calculate length of model paths while scaling flux from the first detector.'
    )

    parser.add_argument('--minscale', type=float, default=102., help='Minimum percentage scale of inst1')
    parser.add_argument('--maxscale', type=float, default=110., help='Maximum percentage scale of inst1.')
    parser.add_argument('-i', '--increment', type=float, default=0.2, help='Scale increments, in percent. Default = 0.2%')
    parser.add_argument('-c', '--choose_order', action='store_true')
    parser.add_argument('--choose_order_scale', type=float, help='Scale used on inst1 for --choose_order run. Default is 0.5*(minscale+maxscale)')
    parser.add_argument('--nolog', action='store_true', help='Do not use log flux for model')
    parser.add_argument('--chi2', action='store_true', default=True, help='In addition to length vs scale, plot X^2 vs scale.')
    parser.add_argument('--models', action='store_true', help='Plot the models generated for each scale factor')
    parser.add_argument('--pdir', default='./data/pickle/scale_lengths', help='Pickle directory. Default is ./data/pickle/scale_lengths')
    parser.add_argument('--pdfdir', default='./plots/scale_lengths', help='PDF directory. Default is ./plots/scale_lengths}')
    parser.add_argument('inst1', help='[lmhs]')
    parser.add_argument('inst2', help='[lmhs]')
    parser.add_argument('epoch', type=int)
    parser.add_argument('rdbs', nargs='+', help='RDB files containing fit results')
    args = parser.parse_args()

    if args.choose_order_scale is None:
        args.choose_order_scale = 0.5*(args.minscale+args.maxscale)

    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()

    inst1 = args.inst1
    inst2 = args.inst2

    # read all data
    hdr, data_rdb = xcal.read_fits_rdbs(args.rdbs)

    # remove epochs without both of the specified instruments
    data_rdb_ = { k : np.array((), dtype=data_rdb[k].dtype) for k in data_rdb }
    groups = xcal.group_days(data_rdb['days'])
    for g in groups:
        ii = np.where((data_rdb['inst'][g]==inst1) |
                      (data_rdb['inst'][g]==inst2))[0]
        if inst1 not in data_rdb['inst'][g]: continue
        if inst2 not in data_rdb['inst'][g]: continue
        for k in data_rdb: data_rdb_[k] = np.append(data_rdb_[k], data_rdb[k][g][ii])
    data_rdb = data_rdb_
    groups = xcal.group_days(data_rdb['days'])

    # retain only specified epoch, if requested
    if args.epoch:
        groups = xcal.group_days(data_rdb['days'])
        if args.epoch > len(groups):
            if rank == 0: sys.stderr.write("there are only {} epochs\n".format(len(groups)))
            exit(1)
        for k in data_rdb: data_rdb[k] = data_rdb[k][groups[args.epoch-1]]

    # for each inst, build arrays of all date/flux data, which are used to build samples
    data = { inst : scale_module.combine_fluxes(data_rdb, inst) for inst in (inst1, inst2) }

    if args.choose_order:
        if rank: return
        choose_order(data, args)
        return

    sim_days = get_sim_days(data, ntimes)

    if rank == 0:
        master(comm, data, sim_days, args, hdr)
    else:
        worker(comm, data, sim_days, args, hdr)

if __name__ == '__main__':
    main()
