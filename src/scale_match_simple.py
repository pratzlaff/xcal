import argparse
import glob
import numpy as np
import sys

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from matplotlib import rc, rcParams
#matplotlib.rc('text', usetex=True)
#rcParams.update({'font.size': 14})

import scale as scale_module
import xcal

# In some cases, all instruments will be combined into one file
# covering a certain wavelength band. In other cases, the instruments
# will be spread across multiple files.
#
# So the data dict will be keyed on wlo
#
def read_rdbs(rdbs):
    d = {}
    for rdbfile in rdbs:
        hdr, data = xcal.read_fits_rdb(rdbfile)
        wlo = hdr['wlo']
        whi = hdr['whi']

        if wlo not in d:
            d[wlo] = {
                'hdr' : hdr,
                'data' : data,
                'wlo' : wlo,
                'whi' : whi,
            }
        else:
            for k in data:
                d[wlo]['data'][k] = np.append(d[wlo]['data'][k], data[k])

    # put times in sequence
    for wlo in d:
        sorti = np.argsort(d[wlo]['data']['days'])
        for k in d[wlo]['data']:
            d[wlo]['data'][k] = d[wlo]['data'][k][sorti]

    hdr = []
    data = []
    wlo = []
    whi = []
    for wlo_ in d:
        hdr.append(d[wlo_]['hdr'])
        data.append(d[wlo_]['data'])
        wlo.append(wlo_)
        whi.append(d[wlo_]['whi'])

    wlo = np.array(wlo).astype(np.float)
    whi = np.array(whi).astype(np.float)
    sorti = np.argsort(wlo)
    data = [data[i] for i in sorti]
    hdr = [hdr[i] for i in sorti]
    wlo = wlo[sorti]
    whi = whi[sorti]
    return hdr, data, wlo, whi

def dist_vs_scale(data, scales, args):
    inst1 = args.inst1
    inst2 = args.inst2

    n1 = data[inst1]['days'].size
    n2 = data[inst2]['days'].size

    days = np.append(data[inst1]['days'], data[inst2]['days'])
    inst = np.append(np.full(n1, inst1, dtype=str), np.full(n2, inst2, dtype=str))
    flux = np.append(data[inst1]['flux'], data[inst2]['flux'])
    err = np.append(data[inst1]['flux_err'], data[inst2]['flux_err'])

    sorti = np.argsort(days)

    days = days[sorti]
    inst = inst[sorti]
    flux = flux[sorti]
    err = err[sorti]

    inst1_mask = inst == inst1

    flux_inst1 = flux[inst1_mask].copy()
    err_inst1 = err[inst1_mask].copy()

    # indices after which the instrument changed
    cpi = np.where((inst[0:-1] != inst[1:]) & (days[1:]-days[0:-1]<0.5))[0]

    dist = np.zeros_like(scales)

    for i in range(scales.size):
        flux[inst1_mask] = scales[i] * flux_inst1
        err[inst1_mask] = scales[i] * err_inst1

        dist[i] = np.sum(np.abs(flux[cpi] - flux[cpi+1]))

    return dist

def main():

    parser = argparse.ArgumentParser(
        description='Simple scale match procedure.'
    )
    parser.add_argument('--minscale', type=float, default=50., help='Minimum percentage scale of inst1')
    parser.add_argument('--maxscale', type=float, default=200., help='Maximum percentage scale of inst1.')
    parser.add_argument('-i', '--increment', type=float, default=0.1, help='Scale increments, in percent. Default = 0.1%')
    parser.add_argument('-p', '--pdf', help='Output PDF file')
    parser.add_argument('inst1', help='[lmhs]')
    parser.add_argument('inst2', help='[lmhs]')
    parser.add_argument('rdbs', nargs='+', help='RDB files containing fit results')
    args = parser.parse_args()

    inst1 = args.inst1
    inst2 = args.inst2

    # read all data
    hdr_all, data_all, wlo_all, whi_all = read_rdbs(args.rdbs)

    scales = 1/100. * (np.arange((args.maxscale-args.minscale)/args.increment+1) * args.increment + args.minscale)

    dist = None

    date_obs = []

    for i in range(len(data_all)):
        hdr = hdr_all[i]
        data = data_all[i]
        wlo = wlo_all[i]
        whi = whi_all[i]

        # remove epochs without both of the specified instruments,
        # retain only fluxes with the two instruments
        data_ = { k : np.array((), dtype=data[k].dtype) for k in data }
        groups = xcal.group_days(data['days'])
        for g in groups:
            ii = np.where((data['inst'][g]==inst1) |
                          (data['inst'][g]==inst2))[0]
            if inst1 not in data['inst'][g]:
                continue
            if inst2 not in data['inst'][g]:
                continue
            for k in data_:
                data_[k] = np.append(data_[k], data[k][g][ii])
        data = data_

        if np.where(data['days'][1:]-data['days'][0:-1]<0)[0].size:
            raise ValueError("times out of sequence")

        groups = xcal.group_days(data['days'])

        data_orig = data

        if dist is None:
            dist = np.zeros((len(data_all), len(groups),  scales.size))

        # each epoch at a time
        for j in range(len(groups)):
            g = groups[j]
            data = { k : data_orig[k][g] for k in data_orig }
            data = { inst : scale_module.combine_fluxes(data, inst) for inst in (inst1, inst2) }
            if i==0:
                date_obs.append(xcal.date_obs_obsid(data[inst1]['obsid'][0]))
            dist[i, j] = dist_vs_scale(data, scales, args)
    #         scale_mins[i, j] = scales[np.argmin(dist[i, j])]

    # np.savetxt(sys.stdout, np.transpose(scale_mins), fmt="%.3f")

    scale_mins = np.zeros((len(data_all) ,len(groups)))

    # for each wavelength band
    for i in range(dist.shape[0]):

        # for each epoch
        for j in range(dist.shape[1]):
            scale_mins[i, j] = scales[np.argmin(dist[i, j])]

    np.savetxt(sys.stdout, np.transpose(scale_mins), fmt="%.3f")

    if args.pdf:
        pdf = PdfPages(args.pdf)
        
    # for each epoch
    for i in range(dist.shape[1]):

        nrows, ncols = 2, 3
        fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(11, 8.5))

        # for each wavelength band
        for j in range(dist.shape[0]):
            row = int(j/ncols) % nrows
            col = j % ncols
            ax = axes[row,col]
            
            plt.sca(axes[row, col])
            ax.autoscale(True)

            x, y = scales, dist[j,i]

            plt.plot(x, y)

            min_at = scale_mins[j,i]
            plt.axvline(x=min_at, color='r', linestyle='-')
            xlim = (min_at-.1, min_at+.1)
            plt.xlim(xlim)

            # FIXME: Y axes aren't being rescaled for the limited X range
            yinterest = y[(x>=xlim[0]) & (x<xlim[1])]
            ylim = (yinterest.min(), yinterest.max())
            plt.ylim(ylim)

            wlo = wlo_all[j]
            whi = whi_all[j]

            ANGSTROM, LAMBDA = "Åλ"
            plt.title('{:.2g}-{:.2g} {}'.format(wlo, whi, ANGSTROM), loc='right')
            if (row==0 and col==0):
                plt.title(date_obs[i][0:10])

            if col==0:
                plt.ylabel('flux NN Abs Dist')

            if row==nrows-1:
                plt.xlabel('HRC-S/LETG Flux Scale')

        plt.tight_layout()
        if args.pdf:
            pdf.savefig()
        else:
            plt.show()
        plt.close()

    if args.pdf:
        pdf.close()
    
if __name__ == '__main__':
    main()
