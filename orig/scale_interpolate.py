import cPickle
import carmcmc as cm
import numpy as np
import argparse

import scale

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from matplotlib import rc
#matplotlib.rc('text', usetex=True)

npaths=1000
nsamples=20000

def get_samples(data, args):

    global nsamples

    models = { }
    samples = { }

    for det in 'leg', 'hrc':

        pickle_dir = 'data/pickle/B%d' % args.band
        pickle_file = '%s/%s_samples.bin' % (pickle_dir, det)

        try:
            f = open(pickle_file, 'r')
            samples[det] = cPickle.load(f)
            f.close()

        except:
            scale.mkdir_p(pickle_dir)
            models[det] = cm.CarmaModel(data[det]['time'], np.log(data[det]['flux']), data[det]['err']/data[det]['flux'], p=4, q=2)
            #    models[det].choose_order(10, njobs=-1) # takes 30 minutes to run on legs
            samples[det] = models[det].run_mcmc(nsamples)
            f = open(pickle_file, 'w')
            cPickle.dump(samples[det], f, 2)
            f.close()

    return samples

def make_plots(data, samples, args):

    global npaths

    band = args.band

    pdf = PdfPages('hrc_scale_interpolate_B%d.pdf' % band)
    fig = plt.figure(figsize = (11, 8.5))
    plot_dims = (2, 3)

    leg_start, leg_stop = scale.get_gaps(data['leg']['time'])
    hrc_start, hrc_stop = scale.get_gaps(data['hrc']['time'])

    for i in range(leg_start.size):

        row = int(i/plot_dims[1])
        col = i % plot_dims[1]
        plt.subplot2grid(plot_dims, (row,col))

        alltimes = np.concatenate((data['leg']['time'][leg_start[i]:leg_stop[i]+1],
                                   data['hrc']['time'][hrc_start[i]:hrc_stop[i]+1]))
        tmax = alltimes.max()
        tmin = alltimes.min()
        span = tmax - tmin

        tmin = tmin - span * .05
        tmax = tmax + span * .05
        span = tmax - tmin

        npoints = 201

        times = np.arange(npoints) * span / npoints + tmin

        sims = get_sims( samples, times, args, plots=True, epoch=i+1)
        mean = { }

        for det in sims:
            sims[det] = np.ma.masked_invalid(sims[det])
            #sims[det] = np.exp(sims[det])
            mean[det] = np.exp(sims[det].mean(axis=0))

        x = data['leg']['time'][leg_start[i]:leg_stop[i]+1]
        y = data['leg']['flux'][leg_start[i]:leg_stop[i]+1]
        yerr = data['leg']['err'][leg_start[i]:leg_stop[i]+1]
        plt.plot(x, y, 'ro', label='ACIS-S')
        plt.plot(times, mean['leg'], 'r--')
        plt.errorbar(x, y, yerr, ecolor='k', fmt=None)

        x = data['hrc']['time'][hrc_start[i]:hrc_stop[i]+1]
        y = data['hrc']['flux'][hrc_start[i]:hrc_stop[i]+1]
        yerr = data['hrc']['err'][hrc_start[i]:hrc_stop[i]+1]
        plt.plot(x, y, 'bs', label='HRC-S')
        plt.plot(times, mean['hrc'], 'b--')
        plt.errorbar(x, y, yerr, ecolor='k', fmt=None)

        plt.xlabel(r'days')
        plt.ylabel('%g - %g keV flux' % scale.band_range(band))

        if i==0: plt.legend()

#        plt.show()

    plt.tight_layout()
    pdf.savefig(fig)
    pdf.close()

def get_sims(samples, times, args, plots=False, epoch=0):

    global npaths

    sims = { }
    sims_shape = ( npaths, times.size )

    for det in samples:

        pickle_dir = 'data/pickle/B%d' % args.band
        if plots:
            pickle_file = '%s/%s_sims_plots_E%d.bin' % (pickle_dir, det, epoch)
        else:
            pickle_file = '%s/%s_sims.bin' % (pickle_dir, det)

        try:
            f = open(pickle_file, 'r')
            sims[det] = cPickle.load(f)
            f.close()

            if cmp( sims[det].shape, sims_shape ):
                raise "The simulations must be recalculated"

        except:
            scale.mkdir_p(pickle_dir)
            sims[det] = np.zeros(sims_shape)
            for i in range(npaths):
                sims[det][i] = samples[det].simulate(times, bestfit='random')

            f = open(pickle_file, 'w')
            cPickle.dump(sims[det], f, 2)
            f.close()

    return sims

def select( data ):

    start = {}
    stop = {}

    indices = np.array([], dtype=np.int)

    for det in data:
        start[det], stop[det] = scale.get_gaps( data[det]['time'] )
    
    for i in range(start['leg'].size):
        a = data['leg']['time'][start['leg'][i]:stop['leg'][i]+1]
        h = data['hrc']['time'][start['hrc'][i]:stop['hrc'][i]+1]

        mask = (h>a.min()) & (h<a.max())

        indices = np.concatenate((indices, np.where(mask)[0] + start['hrc'][i]))

    return indices

def main():

    parser = argparse.ArgumentParser(
        description='Perform CARMA simulation models for ACIS and HRC-S fluxes separately, calculate average of model ratios at HRC-S observation times.'
    )

    parser.add_argument('-b', '--band', type=int, default=0)
    args = parser.parse_args()

    data = { }
    mean = { }

    for det in 'leg', 'hrc':
        rdbfile = 'plots/mkn421_B%d_binned_%s.rdb' % (args.band, det)
        try:
            data[det] = scale.readflux( rdbfile )
        except IOError:
            if rank == 0: sys.stderr.write("could not read file '%s'" % rdbfile)
            exit(1)

    indices = select( data )

#    for k,v in data['hrc']: v = v[indices]
    for key in data['hrc']: data['hrc'][key] = data['hrc'][key][indices]

    samples = get_samples( data, args )

    sims = get_sims( samples, data['hrc']['time'], args )

    for det in sims:
        mean[det] = np.exp(sims[det].mean(axis=0))

    ratio = mean['leg'] / mean['hrc']
    ratio_mean = ratio.mean()
    ratio_std = ratio.std() / np.sqrt(npaths)

    print("ratio mean = %f (%f)" % ( ratio_mean, ratio_std ) )

    make_plots(data, samples, args)

    pass

if __name__ == '__main__':
    main()
