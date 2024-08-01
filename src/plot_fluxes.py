import pprint
import argparse
import re
import sys

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from matplotlib import rc, rcParams
rc('text', usetex=True)
rcParams.update({'font.size': 14})

import xcal

titles = { 'mkn421' : 'Mkn 421', 'pks2155' : 'PKS 2155' }

def read_rdb(rdbfile):
    params = {}
    try:
        fh = open(rdbfile)
        for line in fh:
            match = re.match(r'^#\s+(\w+):\s+(\S+)', line)
            if match:
                params[match.groups()[0]] = match.groups()[1]
            else:
                jnk = fh.readline()
                obsid, binn, days, inst, order, flux, err = np.loadtxt(fh, dtype='int,int,float,U1,U3,float,float', unpack=True)
    finally:
        fh.close()

    flux, err = np.ma.masked_array(flux), np.ma.masked_array(err)

    return params, obsid, binn, days, inst, order, flux, err

def nan_filter(flux, err):
    return np.where((flux==flux) & (err==err))[0]

def median_filter(err):
    return np.where(err < 20*np.median(err))[0]

def filter_data(args, inst):
    ind = np.array((), dtype=int)
    for c in args.inst:
        ind = np.append(ind, np.where(c == inst))
    sorti = np.argsort(ind)
    return ind[sorti]

def group_days(days):
    groups = []
    this_group=[0]
    # if more than three days past the last entry, start over
    for i in range(1, days.size):
        if (days[i]-days[i-1])>3:
            groups.append(this_group)
            this_group = [i]
        else:
            this_group.append(i)
    groups.append(this_group)
    return groups

def construct_group(args, obsid, days, inst, order, flux, err):
    # return a dict of the form
    # { 'l' : { 'neg' : { obsid : [],
    #                     days  : [],
    #                     flux : [],
    #                     err : []
    #                   }
    #           'pos' : { ... }
    #         }
    # }
    group = {}
    for c in args.inst:
        for o in ('neg', 'pos'):
            ind = np.where((c == inst) & (o==order))[0]
            if ind.size:
                if not c in group: group[c] = {}
                group[c][o] = {
                    'obsid' : obsid[ind],
                    'days' : days[ind],
                    'flux' : flux[ind],
                    'err' : err[ind]
                }
    return group

def get_groups(args):

    params, obsid, binn, days, inst, order, flux, err = read_rdb(args.rdb)
    ind = filter_data(args, inst)

    obsid = obsid[ind]
    binn = binn[ind]
    days = days[ind]
    inst = inst[ind]
    order = order[ind]
    flux = flux[ind]
    err = err[ind]

    n = err.size
    ind = nan_filter(flux, err)
    sys.stderr.write("NaN filter retained {} / {} values\n".format(ind.size, n))
    obsid = obsid[ind]
    binn = binn[ind]
    days = days[ind]
    inst = inst[ind]
    order = order[ind]
    flux = flux[ind]
    err = err[ind]

    if args.filter:
        n = err.size
        ind = median_filter(err)
        sys.stderr.write("Median error filter retained {} / {} values\n".format(ind.size, n))
        obsid = obsid[ind]
        binn = binn[ind]
        days = days[ind]
        inst = inst[ind]
        order = order[ind]
        flux = flux[ind]
        err = err[ind]

    # each group is a set of interleaved observations
    group_ind = xcal.group_days(days)

    groups = [ ]

    for gi in group_ind:
        group = construct_group(args,
                                obsid[gi],
                                days[gi],
                                inst[gi],
                                order[gi],
                                flux[gi],
                                err[gi]
        )
                      
        if len(group.keys())==1 and not args.all:
            continue
        groups.append(group)

    #pprint.pprint(groups)
    return groups, params

def plot(args):

    groups, params = get_groups(args)

    if args.pdf:
        pdf = PdfPages(args.pdf)
        fig = plt.figure(figsize = (11, 8.5))

    legend=True

    plot_dims = (1, 1)

    for i in range(len(groups)):

        row = int(i/plot_dims[1]) % plot_dims[0]
        col = i % plot_dims[1]
        plt.subplot2grid(plot_dims, (row,col))

        if not row and not col:
            instruments = list(groups[i].keys())
            orders = list(groups[i][instruments[0]].keys())
            date_obs = xcal.date_obs_obsid(groups[i][instruments[0]][orders[0]]['obsid'][0])
            plt.title(r'\textrm{' + titles[params['source']] + ': ' + date_obs[:10] +'}')

        # if not row and col==1:
        #     plt.title(r'\textrm{CalDB 4.7.4}')

        if row == plot_dims[0]-1:
            plt.xlabel(r'\textrm{Time Offset (ks)}')

        if col == 0:
            elow, ehigh = 12.39854/float(params['wavmax']), 12.39854/float(params['wavmin'])
            ylabel = r'$$\textrm{'+'{:.2g}-{:.2g} keV Flux '.format(elow,ehigh)+r'(erg s}^{-1}\textrm{ cm}^{-2}\textrm{)}$$'
            plt.ylabel(ylabel)

        plot_group(groups[i], args, params, legend=legend)
        if i == 0:
            plt.legend()

        legend=False

        if (
                (row == plot_dims[0]-1 and col==plot_dims[1]-1) or
                (i==len(groups)-1)
        ):
            plt.tight_layout()
            if args.pdf:
                pdf.savefig(fig)
            else:
                plt.show()
            plt.clf()

    if args.pdf:
        pdf.close()

def group_dminmax(group):
    dmin, dmax = None, None
    for inst in group:
        for order in group[inst]:
            dmin_ = group[inst][order]['days'][0]
            dmax_ = group[inst][order]['days'][-1]

            if dmin is None:
                dmin = dmin_
            else:
                if dmin_ < dmin:
                    dmin=dmin_

            if dmax is None:
                dmax = dmax_
            else:
                if dmax_ > dmax:
                    dmax=dmax_

    return dmin, dmax

def group_days_offset(group):
    dmin, dmax = group_dminmax(group)
    return dmin

def plot_group(group, args, params, legend=False):
    global titles
    
    # FIXME: offsets are from the middle of the first bin, rather than
    # from the start of the first bin
    dmin, dmax = group_dminmax(group)

    colors = {'s' : 'b', 'l' : 'k', 'm' : 'r', 'h' : 'm' }
    labels = {
        's' : 'HRC/LEG',
        'l' : 'ACIS/LEG',
        'm' : 'ACIS/MEG',
        'h' : 'ACIS/HEG',
    }

    markers = { 'neg' : 's' , 'pos' : 'o' }
    linestyles = { 'neg' : '--' , 'pos' : '-' } # for the error bars

    tfudge_mult = { 'neg' : -1, 'pos' : +1 }
    for inst in group:
        for order in group[inst]:
            time = 86400 * (group[inst][order]['days'] - dmin) / 1e3
            tfudge = tfudge_mult[order] * 0.003 * (time.max() - time.min())
            time += tfudge
            #plt.plot(time, group[inst][order]['flux'], colors[inst]+markers[order], linestyle='None')
            color = colors[inst]
            marker = markers[order]
            label = None
            if order=='pos':
                label = labels[inst]
            eb = plt.errorbar(time,
                              group[inst][order]['flux'],
                              yerr=group[inst][order]['err'],
                              fmt=colors[inst]+'.',#markers[order],
                              label=label,
            )
            # eb[-1][0] is the LineCollection objects of the errorbar lines
            eb[-1][0].set_linestyle(linestyles[order])

                
def main():

    parser = argparse.ArgumentParser(
        description='Plot observed fluxes from interleaved calibration observations.'
    )
    parser.add_argument('-p', '--pdf', help='Save plot to named file.')
    parser.add_argument('-f', '--filter', action='store_true', help='Median error filter.')
    parser.add_argument('-a', '--all', action='store_true', help='Also plot series which only have one instrument configuration.')
    parser.add_argument('rdb', help='RDB file')
    parser.add_argument('inst', help='Instrument code (e.g., "sl" for HRC-S/LEG and ACIS-S/LEG). s=HRC-S/LEG, l=ACIS-S/LEG, m=ACIS-S/MEG, h=ACIS-S/HEG')
    args = parser.parse_args()

    plot(args)

if __name__ == '__main__':
    main()
