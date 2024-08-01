import argparse
import carmcmc as cm
import cPickle
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from matplotlib import rc, rcParams
#rc('text', usetex=True)
#rcParams.update({'font.size': 14})

import scale
import xcal

npaths=1000
nsamples=20000

def get_samples(data, args, hdr):

    global nsamples

    models = { }
    samples = { }

    for inst in 'l', 's':

        pickle_file = '{}/{}_{}_B{}_samples.bin'.format(args.pdir, hdr['source'], inst, hdr['band'])

        try:
            f = open(pickle_file, 'r')
            samples[inst] = cPickle.load(f)
            f.close()

        except:
            xcal.mkdir_p(args.pdir)
            models[inst] = cm.CarmaModel(data[inst]['days'], np.log(data[inst]['flux']), data[inst]['err']/data[inst]['flux'], p=4, q=2)
            #    models[inst].choose_order(10, njobs=-1) # takes 30 minutes to run on legs
            samples[inst] = models[inst].run_mcmc(nsamples)
            f = open(pickle_file, 'w')
            cPickle.dump(samples[inst], f, 2)
            f.close()

    return samples

def make_plots(data, samples, args, hdr):

    global npaths

    group_l = xcal.group_days(data['l']['days'])
    group_s = xcal.group_days(data['s']['days'])

    if args.pdf:
        pdf = PdfPages(args.pdf)
        fig = plt.figure(figsize = (8.5, 11))

    if len(group_l) < 5:
        plot_dims = (len(group_l), 1)
    else:
        plot_dims = (4,2)

    if len(group_l) != len(group_s):
        raise ValueError('ACIS/LETG and HRC/LETG do not have the same number of groups')

    for i in range(len(group_l)):
        gl = group_l[i]
        gs = group_s[i]

        row = int(i/plot_dims[1]) % plot_dims[0]
        col = i % plot_dims[1]
        plt.subplot2grid(plot_dims, (row,col))

        date_obs = xcal.date_obs_obsid(data['l']['obsid'][gl][0])[:10]

        alldays = np.append(data['l']['days'][gl], data['s']['days'][gs])
        dmax = alldays.max()
        dmin = alldays.min()
        span = dmax - dmin

        dmin = dmin - span * .05
        dmax = dmax + span * .05
        span = dmax - dmin

        npoints = 201

        days = np.arange(npoints) * span / npoints + dmin

        sims = get_sims(samples, days, args, hdr, plots=True, epoch=i+1)
        mean = { inst : np.exp(sims[inst].mean(axis=0)) for inst in sims }

        title = r'\textrm{' + date_obs[:10] +'}'
        title = date_obs[:10]
        if not cmp((row,col),(0,0)):
            titles = { 'mkn421' : 'Mkn 421', 'pks2155' : 'PKS 2155' }
            title = r'\textrm{' + titles[hdr['source']] + ': ' + date_obs[:10] + '}'
            title = titles[hdr['source']] + ': ' + date_obs[:10]
            
        plt.title(title)

        x = data['l']['days'][gl]-dmin
        y = data['l']['flux'][gl]
        yerr = data['l']['err'][gl]
        plt.plot(x, y, 'ro', label='ACIS-S')
        plt.plot(days-dmin, mean['l'], 'r--')
        plt.errorbar(x, y, yerr, ecolor='k', fmt='none')

        x = data['s']['days'][gs]-dmin
        y = data['s']['flux'][gs]
        yerr = data['s']['err'][gs]
        plt.plot(x, y, 'bs', label='HRC-S')
        plt.plot(days-dmin, mean['s'], 'b--')
        plt.errorbar(x, y, yerr, ecolor='k', fmt='none')

        if row==plot_dims[0]-1:
            plt.xlabel(r'\textrm{Days}')
            plt.xlabel('Days')

        elo, ehi = xcal.band_energy_limits(int(hdr['band']))
        ylabel = r'$$\textrm{'+'{:.2g}-{:.2g} keV Flux '.format(elo,ehi) + r'}$$'
        ylabel = '{:.2g}-{:.2g} keV Flux '.format(elo,ehi)
        plt.ylabel(ylabel)

        if i==0: plt.legend()

        if not cmp((row+1,col+1), plot_dims) or i==len(group_l)-1:
            plt.tight_layout()
            if args.pdf:
                pdf.savefig(fig)
                plt.clf()
            else:
                plt.show()

    if args.pdf:
        pdf.close()

def get_sims(samples, days, args, hdr, plots=False, epoch=0):

    global npaths

    sims = { }
    sims_shape = ( npaths, days.size )

    for inst in samples:

        if plots:
            pickle_file = '{}/{}_{}_B{}_E{:02d}_sims_plots.bin'.format(args.pdir, hdr['source'], inst, hdr['band'], epoch)
        else:
            pickle_file = '{}/{}_{}_B{}_sims.bin'.format(args.pdir, hdr['source'], inst, hdr['band'])

        try:
            f = open(pickle_file, 'r')
            sims[inst] = cPickle.load(f)
            f.close()

            if cmp( sims[inst].shape, sims_shape ):
                raise "The simulations must be recalculated"

        except:
            xcal.mkdir_p(args.pdir)
            sims[inst] = np.zeros(sims_shape)
            for i in range(npaths):
                sims[inst][i] = samples[inst].simulate(days, bestfit='random')

            f = open(pickle_file, 'w')
            cPickle.dump(sims[inst], f, 2)
            f.close()

    sims = { k : np.ma.masked_invalid(sims[k]) for k in sims }
    return sims

def scale_interpolate(args):

    # all data...
    hdr, data = xcal.read_fits_rdbs(args.rdbs)

    # filter out HRC data which aren't sandwiched between ACIS observations
    data_ = { k : np.array((), dtype=data[k].dtype) for k in data }
    groups = xcal.group_days(data['days'])
    for g in groups:
        ii = np.where(data['inst'][g]=='l')[0]
        if not ii.size: continue
        dmin = data['days'][g][ii].min()
        dmax = data['days'][g][ii].max()
        ii = np.where((data['inst'][g]!='s') |
                              (
                                  (data['days'][g]>dmin) & (data['days'][g]<dmax)
                              ))[0]
        for k in data: data_[k] = np.append(data_[k], data[k][g][ii])
    data = data_

    # filtering out HETG, use only those groups with both ACIS and HRC
    # observations
    data_ = { k : np.array((), dtype=data[k].dtype) for k in data }
    groups = xcal.group_days(data['days'])
    for g in groups:
        ii = np.where((data['inst'][g]=='l')|(data['inst'][g]=='s'))[0]
        if 'l' not in data['inst'][g]: continue
        if 's' not in data['inst'][g]: continue
        for k in data: data_[k] = np.append(data_[k], data[k][g][ii])
    data = data_

    # for each detector, build arrays of all date/flux data, which are used to build samples
    data_samp = { inst : scale.combine_fluxes(data, inst) for inst in ('l', 's') }
    samples = get_samples(data_samp, args, hdr)

    sims = get_sims(samples, data_samp['s']['days'], args, hdr)

    mean = { inst : np.exp(sims[inst].mean(axis=0)) for inst in data_samp }

    ratio = mean['l'] / mean['s']
    ratio_mean = ratio.mean()
    ratio_std = ratio.std() / np.sqrt(npaths)

    print("ratio mean = %f (%f)" % ( ratio_mean, ratio_std ) )

    make_plots(data_samp, samples, args, hdr)

def main():
    parser = argparse.ArgumentParser(
        description='Perform CARMA simulation models for ACIS and HRC LETG fluxes separately, calculate average of model ratios at HRC observation times.'
    )
    parser.add_argument('-p', '--pdf', help='Save plot to named file.')
    parser.add_argument('--pdir', default='./data/pickle/scale_interpolate', help='Pickle directory. Default is ./data/pickle/scale_interpolate')
    parser.add_argument('rdbs', nargs='+', help='RDB files containing fit results')
    args = parser.parse_args()
    scale_interpolate(args)

if __name__ == '__main__':
    main()
