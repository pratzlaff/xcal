import cPickle
import carmcmc as cm
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import matplotlib
#matplotlib.rc('text', usetex=True)

npaths=1000
nsamples=20000
ntimes=201

def get_sim_times(time):

    global ntimes

    sorti = np.argsort(time)
    time = time[sorti]

    start, stop = get_gaps(time)

    sim_times = np.zeros((start.size, ntimes))

    for i in range(start.size):
        time_i = time[start[i]:stop[i]+1]
        tmin = time_i.min()
        tmax = time_i.max()
        span = tmax-tmin
        sim_times[i]= np.arange(ntimes) * span / ntimes + tmin

    return sim_times

def get_sims(sample, times, scale):

    global npaths

    pickle_file = 'data/pickle/combined_sims_%04.1f.bin' % (100*(scale-1))

    sims_shape = ( times.shape[0], npaths, times.shape[1] )
    try:
        f = open(pickle_file, 'r')
        sims = cPickle.load(f)
        f.close()

        if cmp(sims_shape, sims.shape):
            raise "The simulations must be recalculated"

    except:
        sims = np.zeros(sims_shape)
        for i in range(sims.shape[0]):
            for j in range(sims.shape[1]):
                sims[i,j] = sample.simulate(times[i], bestfit='random')

        f = open(pickle_file, 'w')
        cPickle.dump(sims, f, 2)
        f.close()

    return sims

def get_samples(time, flux, err, scale):

    global nsamples

    pickle_file = 'data/pickle/combined_samples_%04.1f.bin' % (100*(scale-1))

    try:
        f = open(pickle_file, 'r')
        sample = cPickle.load(f)
        f.close()

    except:
        model = cm.CarmaModel(time, np.log(flux), err/flux, p=7, q=5)
        #    model.choose_order(10, njobs=-1) # takes 30 minutes to run on legs
        sample = model.run_mcmc(nsamples)
        f = open(pickle_file, 'w')
        cPickle.dump(sample, f, 2)
        f.close()

    return sample

def get_gaps(time):
    gaps, = np.where(time[1:] - time[0:-1] > 10)
    splits = np.split(time, gaps+1)

    start = np.zeros(gaps.size+1, dtype=np.long)
    stop = start.copy()

    start[1:] = gaps+1
    stop[0:-1] = gaps ; stop[-1] = time.size-1

    return start, stop

def readflux(file):
    obsid, days, pos, pos_low, pos_high, neg, neg_low, neg_high = np.loadtxt(file,
        skiprows=10, unpack=True )
    flux = pos + neg
    flux_err = 0.5 * np.sqrt( (pos_high-pos_low)**2 + (neg_high-neg_low)**2 )

    return { 'time':days, 'flux':flux, 'err':flux_err }

if __name__ == '__main__':

    corr = np.arange(71) * 10 / 100. + 2.
#    corr = np.arange(5) + 5.0
    scale = 1 + corr/100.

    acis = readflux( 'plots/mkn421_B0_binned_acis_letg.rdb' )
    hrc = readflux( 'plots/mkn421_B0_binned_hrc_letg.rdb' )

    acis_start, acis_stop = get_gaps(acis['time'])
    hrc_start, hrc_stop = get_gaps(hrc['time'])

    time = np.concatenate((acis['time'], hrc['time']))
    sim_times = get_sim_times(time)

    dist = np.zeros(scale.shape)

    pdf = PdfPages('plots_scaled.pdf')
    fig = plt.figure(figsize = (11, 8.5))

    for i in range(scale.size):

        flux = np.concatenate((acis['flux'], hrc['flux'] * scale[i]))
        err = np.concatenate((acis['err'], hrc['err'] * scale[i]))

        sample = get_samples(time, flux, err, scale[i])

        sims = get_sims(sample, sim_times, scale[i])

        # a few scattered NaN are popping up in the simulations, and
        # they screw up everything afterwards
        sims = np.ma.masked_invalid(sims)

        # makes little difference whether we exponentiate before or after taking the mean
        mean = np.exp(sims.mean(axis=1))

        # makes little difference which of these are used
        dist[i] = np.sqrt((mean[:,1:] - mean[:,:-1])**2).sum()
        dist[i] = np.sqrt((1e8*(mean[:,1:] - mean[:,:-1]))**2 + (sim_times[:,1:] - sim_times[:,:-1])**2).sum()


        plot_dims = (2, 2)
        for j in range(acis_start.size):
            row = int( j / plot_dims[1] )
            col = j % plot_dims[1]
            plt.subplot2grid(plot_dims, (row,col))

            x = acis['time'][acis_start[j]:acis_stop[j]+1] - sim_times[j].min()
            y = acis['flux'][acis_start[j]:acis_stop[j]+1]
            yerr = acis['err'][acis_start[j]:acis_stop[j]+1]
            plt.plot(x, y, 'r.', label='ACIS-S')
            plt.errorbar(x, y, yerr, ecolor='k', fmt=None)

            x = hrc['time'][hrc_start[j]:hrc_stop[j]+1] - sim_times[j].min()
            y = hrc['flux'][hrc_start[j]:hrc_stop[j]+1] * scale[i]
            yerr = hrc['err'][hrc_start[j]:hrc_stop[j]+1] * scale[i]
            plt.plot(x, y, 'b.', label='scaled HRC-S')
            plt.errorbar(x, y, yerr, ecolor='k', fmt=None)

#            plt.xlabel(r'Days since 1998.0')
            if not j%2: plt.ylabel(r'Flux (ergs / s / cm^2)')

            plt.plot(sim_times[j]-sim_times[j].min(), mean[j], 'k-')

            if j==0:
                plt.title("scale = %.3f" % scale[i])
                plt.legend(loc='lower left', fontsize=8)

        plt.tight_layout()
        pdf.savefig(fig)

        print "%g\t%g" % (scale[i], dist[i])

        pass

    if (j != np.prod(plot_dims)-1):
        plt.tight_layout()
        pdf.savefig(fig)
    pdf.close()

    plt.close()

    pdf = PdfPages('lengths.pdf')
    fig = plt.figure(figsize = (11, 8.5))
    plt.plot(scale, dist)
    plt.xlabel('HRC-S scale factor')
    plt.ylabel('combined model path length')
    plt.tight_layout()
    pdf.savefig(fig)
    pdf.close()

    plt.close()

    # plt.plot(scale, dist)
    # plt.xlabel('HRC-S scale factor')
    # plt.ylabel('model distance')
    # plt.show()

