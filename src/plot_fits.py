import argparse
import pprint
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import re

import xcal

from matplotlib import rc, rcParams
rc('text', usetex=True)
rcParams.update({'font.size': 14})

titles = { 'mkn421' : 'Mkn 421', 'pks2155' : 'PKS 2155' }

group_number=0

def plot(args):

    hdr, data = xcal.read_fits_rdbs(args.rdbs)

    groups = xcal.group_days(data['days'])

    if args.pdf:
        pdf = PdfPages(args.pdf)
        dims = (8.5, 11)
        try:
            if int(hdr['band']) == 0:
                dims = (11, 8.5)
        except:
            pass
        fig = plt.figure(figsize = dims)

    legend=True
    for i in range(len(groups)):

        if np.unique(data['inst'][groups[i]]).size == 1 and not args.all:
            continue

        plot_group(hdr, data, groups[i], args, legend=legend)
        legend=False

        plt.tight_layout()
        if args.pdf:
            pdf.savefig(fig)
            plt.clf()
        else:
            plt.show()

    if args.pdf:
        pdf.close()

def print_stds(hdr, data, args):
    global group_number

    if not args.sigma:
        return

    if not group_number:
        print("# group\tinst\tparam\tsigma")
        #print("# group\tinst\torders\tflux\tflux_low\tflux_high\t")

    group_number += 1

    for param in f['HRC']['LEG']['neg']:
        for det in f:
            for part in f[det]:
                inst = det + '/' + part

                best_all = list(f[det][part]['neg'][param]['best']) # ensure we get a copy
                best_all.extend(f[det][part]['pos'][param]['best'])
                print('{}\t{}\t{}\t{:g}'.format(group_number, param, inst, np.std(best_all, ddof=1)))


            # for orders in 'neg', 'pos':
            #     flux = np.array(f[det][part][orders]['flux']['best'])
            #     flux_low = np.array(f[det][part][orders]['flux']['min'])
            #     flux_high = np.array(f[det][part][orders]['flux']['max'])
            #     for i in range(flux.size):
            #         print('{}\t{}\t{}\t{:g}\t{:g}\t{:g}'.format(group_number, inst, orders, flux[i], flux_low[i], flux_high[i]))
            

def plot_group(hdr, data, group, args, legend=False):
    global titles

    print_stds(hdr, data, args)
    
    doff = data['days'][group][0]
    days_ = data['days'][group] - doff
    dfudge = .003 * days_[-1]

    colors = {
        'l' : 'k',
        'm' : 'r',
        'h' : 'm',
        's' : 'b'
    }

    markers = { 'neg' : 'v' , 'pos' : '^' }
    linestyles = { 'neg' : '--' , 'pos' : '-' } # for the error bars

    params = hdr['params']

    try:
        wlo, whi = xcal.band_wav_limits(int(hdr['band']))
    except:
        wlo, whi = hdr['wlo'], hdr['whi']
        
    labels={
        'alpha' : r'$$\alpha$$',
        'beta' : r'$$\beta$$',
        'PhoIndex' : r'$$\alpha$$',
        'norm' : r'\textrm{Norm}',
        'flux' : r'$$\textrm{'+'{:.2g}-{:.2g} \AA Flux '.format(float(wlo), float(whi))+r'(erg s}^{-1}\textrm{ cm}^{-2}\textrm{)}$$'
        }

    if len(params) >= 4:
        plot_dims = (2, int((len(params)+1)/2))
    else:
        plot_dims = (len(params), 1)

    date_obs = xcal.date_obs_obsid(data['obsid'][group][0])

    for i in range(len(params)):
        param = params[i]
        row = int(i/plot_dims[1]) % plot_dims[0]
        col = i % plot_dims[1]
        plt.subplot2grid(plot_dims, (row,col))

        if not row and not col:
            plt.title(r'\textrm{' + titles[hdr['source']] + ': ' + date_obs[:10] +'}')

        if not row and col==1:
            plt.title(r'\textrm{CalDB 4.7.4}')

        if row == plot_dims[0]-1:
            plt.xlabel(r'\textrm{Time Offset (ks)}')

        days = { 'neg' : days_ - dfudge, 'pos' : days_ + dfudge }
        for inst in colors:

            for orders in 'neg', 'pos':
                ii = np.where(
                    (data['inst'][group]==inst) & (data['order'][group]==orders)
                )[0]

                if not ii.size:
                    continue

                best = data[param][group][ii]
                low = data[param+'_min'][group][ii]
                high = data[param+'_max'][group][ii]
                                    
                label=None
                if orders=='pos':
                    det = list(xcal.configs[inst].keys())[0]
                    part = xcal.configs[inst][det]
                    label=r'\textrm{' + det+'/'+part + '}'
                        
                yerr=None
                if not args.noerrb:
                    yerr = [best-low, high-best]
                    yerravg = 0.5*(yerr[1]+yerr[0])
                    mask = yerravg < 5*np.median(yerravg) 
                    eb = plt.errorbar(days[orders][ii][mask]*86.4,
                                      best[mask],
                                      yerr=(yerr[i][mask] for i in range(len(yerr))),
                                      ecolor=colors[inst],
                                      fmt=colors[inst]+markers[orders],
                                      capthick=0,
                                      label=label
                    )

                    # eb[-1][0] is the LineCollection objects of the errorbar lines
                    eb[-1][0].set_linestyle(linestyles[orders])

                else:
                    plt.plot(days[orders][ii]*86.4, best, colors[inst]+markers[orders], label=label)

                plt.ylabel(labels.get(param, param))

        if legend and not row and not col:
            plt.legend(loc='upper left', fontsize=10)
                
def main():

    parser = argparse.ArgumentParser(
        description='Plot interleaved calibration observation fit results.'
    )
    parser.add_argument('-p', '--pdf', help='Save plot to named file.')
    parser.add_argument('--noerrb', action='store_true', help='Do not plot error bars.')
    parser.add_argument('-s', '--sigma', action='store_true', help='Print np.std(ddof=1) of each param best-fit value.')
    parser.add_argument('-a', '--all', action='store_true', help='Plot groups with only one detector.')
    parser.add_argument('rdbs', nargs='+', help='RDB files containing fit results')
    args = parser.parse_args()

    plot(args)

if __name__ == '__main__':
    main()
