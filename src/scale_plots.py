import argparse
import glob
import math
import numpy as np
import sys

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from matplotlib import rc, rcParams
#matplotlib.rc('text', usetex=True)
#rcParams.update({'font.size': 14})

import scale as scale_module
import xcal

def plot(data, scales, args, hdr):

    epoch = args.epoch
    source = hdr['source']
    i1 = args.inst1
    i2 = args.inst2

    groups = { inst : xcal.group_days(data[inst]['days']) for inst in data }

    if args.pdf:
        pdf = PdfPages(args.pdf)
        fig = plt.figure(figsize = (11, 8.5))

    ngroups = len(groups[list(groups.keys())[0]])

    plot_dims = (1, 1)

    labels = { 'l': 'ACIS/LEG', 
               'm': 'ACIS/MEG',
               'h': 'ACIS/HEG',
               's': 'HRC/LEG' }

    symbols = { 'l': 'ko',
                'm': 'ro',
                'h': 'mo',
                's': 'bo' }

    colors = { 'l': 'k',
               'm': 'r',
               'h': 'm',
               's': 'b' }

    for i in range(scales.size):

        for j in range(ngroups):

            date_min = None
            for inst in data:
                g = groups[inst][j]
                min = data[inst]['days'][g].min()
                if date_min is None or min < date_min: date_min = min

            row = int( j / plot_dims[1] )
            col = j % plot_dims[1]

            plt.subplot2grid(plot_dims, (row,col))

            for inst in data:
                g = groups[inst][j]

                x = data[inst]['days'][g] - date_min
                y = data[inst]['flux'][g]
                y_err = data[inst]['flux_err'][g]
                y_neg = data[inst]['flux_neg'][g]
                y_neg_err = data[inst]['flux_neg_err'][g]
                y_pos = data[inst]['flux_pos'][g]
                y_pos_err = data[inst]['flux_pos_err'][g]

                label = labels[inst]

                if inst == i1:
                    y = y * scales[i]
                    y_err = y_err * scales[i]
                    y_neg = y_neg * scales[i]
                    y_neg_err = y_neg_err * scales[i]
                    y_pos = y_pos * scales[i]
                    y_pos_err = y_pos_err * scales[i]
                    label = 'Scaled ' + label

                if args.noerrbars:
                    y_err = None
                    y_pos_err = None
                    y_neg_err = None

                #ind1 = np.where(y_err < 5*np.median(y_err))[0]
                #ind2 = np.where(y_err >= 5*np.median(y_err))[0]
                #plt.errorbar(x[ind1], y[ind1], y_err[ind1], fmt=symbols[inst], label=label)
                #plt.plot(x[ind2], y[ind2], symbols[inst])

                # mask = y_err < 5*np.median(y_err)
                # plt.errorbar(x[mask], y[mask], y_err[mask], fmt=symbols[inst], label=label)
                # plt.plot(x[~mask], y[~mask], symbols[inst])

                plt.plot(x, y_neg, colors[inst]+'v', label=label)
                plt.plot(x, y_pos, colors[inst]+'^')


            if row==plot_dims[0]-1:
                plt.xlabel(r'Days')
            if col==0:
                try:
                    wlo, whi = xcal.band_wav_limits(int(hdr['band']))
                except:
                    try:
                        wlo, whi = hdr['wlo'], hdr['whi']
                    except:
                        wlo, whi = hdr['wavmin'], hdr['wavmax']
                plt.ylabel('{:g} - {:g}'.format(float(wlo), float(whi)) + r' $\AA$ Flux')

            if j==0:
                plt.title("Scale = %.3f" % scales[i])
                plt.legend(loc='lower left', fontsize=8)

        plt.tight_layout()
        if args.pdf:
            pdf.savefig(fig)
        else:
            plt.show()
        plt.clf()

    if args.pdf:
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

    parser = argparse.ArgumentParser(
        description='Perform CARMA simulation models for two detectors\' combined fluxes, calculate length of model paths while scaling flux from the first detector.'
    )

    parser.add_argument('--minscale', type=float, default=102., help='Minimum percentage scale of inst1')
    parser.add_argument('-p', '--pdf', help='Output PDF file')
    parser.add_argument('--maxscale', type=float, default=110., help='Maximum percentage scale of inst1.')
    parser.add_argument('-i', '--increment', type=float, default=0.2, help='Scale increments, in percent. Default = 0.2%')
    parser.add_argument('-n', '--noerrbars', action='store_true', help='Do not plot error bars.')
    parser.add_argument('inst1', help='[lmhs]')
    parser.add_argument('inst2', help='[lmhs]')
    parser.add_argument('epoch', type=int)
    parser.add_argument('rdbs', nargs='+', help='RDB files containing fit results')
    args = parser.parse_args()

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

    # retain only specified epoch
    groups = xcal.group_days(data_rdb['days'])
    if args.epoch > len(groups):
        sys.stderr.write("there are only {} epochs\n".format(len(groups)))
        exit(1)
    for k in data_rdb: data_rdb[k] = data_rdb[k][groups[args.epoch-1]]

    # for each inst, build arrays of all date/flux data, which are used to build samples
    data = { inst : scale_module.combine_fluxes(data_rdb, inst) for inst in (inst1, inst2) }

    scales = 1/100. * (np.arange((args.maxscale-args.minscale)/args.increment+1) * args.increment + args.minscale)
    plot(data, scales, args, hdr)

if __name__ == '__main__':
    main()
