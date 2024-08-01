import numpy as np
import argparse
import sys

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import rc, rcParams
rc('text', usetex=True)
#rcParams.update({'font.size': 14})

import response
import xcal
import util
import flux

header_printed = False

def print_header(fit, args):
    print('# source: {}'.format(args.source))
    if args.band is not None:
        print('# band: {}'.format(args.band))
    else:
        print('# wlo: {}'.format(args.wlo))
        print('# whi: {}'.format(args.whi))
    cols = ['obsid', 'bin', 'days', 'inst', 'order']
    types = ['N', 'N', 'N', 'S', 'S']

    params = []
    for key in fit:
        params.extend((key, key+'_min', key+'_max'))

    cols.extend(params)
    types.extend(['N']*len(params))

    print('\t'.join(cols))
    print('\t'.join(types))

def print_fit(fit, args, obsid, bin, days, part, order):

    global header_printed
    if not header_printed:
        print_header(fit, args)
        header_printed = True

    code = xcal.configs[args.det.upper()][part.upper()]

    cols = ['{}'.format(d) for d in (obsid, bin)]
    cols.extend(['{:.5f}'.format(days)])
    cols.extend([code, order])

    for param in fit:
        cols.extend(['{:.4g}'.format(fit[param][s]) for s in ('best', 'min', 'max')])
    print('\t'.join(cols))

def handle_obsid(obsid, args):

    global header_printed

    detnam, grating = xcal.instruments_obsid(obsid)
    det = { 'ACIS' : 'acis', 'HRC-' : 'hrc' }[detnam[0:4]]
    if det != args.det:
        return

    tmin, tmax, tmin_binned, tmax_binned = xcal.trange_obsid(obsid)

    bin = [i for i in range(-1, len(tmin_binned))]
    tmin = [tmin]
    tmin.extend(tmin_binned)
    tmax = [tmax]
    tmax.extend(tmax_binned)

    orders = ['neg', 'pos']
    parts = { 'LETG' : ('LEG',), 'HETG' : ('MEG', 'HEG') }[grating]

    # special cases
    if args.band == 1 and grating == 'HETG':
        return
    if args.band == 2 and grating == 'HETG':
        parts = ('MEG',)

    fits = { }
    for part in parts:
        fits[part] = {}
        for order in orders:
            try:
                fit_main, fit_bin = xcal.get_fits(obsid, part, order, args.band, args.wlo, args.whi)
                fits_ = [fit_main]
                fits_.extend(fit_bin)
                fits[part][order] = fits_
            except:
                if args.band is not None:
                    sys.stderr.write("Failed to get fits for {}_{}_{}_B{}\n".format(obsid, part, order, args.band))
                else:
                    sys.stderr.write("Failed to get fits for {}_{}_{}_w{:g}-{:g}\n".format(obsid, part, order, args.wlo, args.whi))

    for i in range(len(tmin)):
        for part in parts:
            for order in orders:
                days = 0.5*(tmax[i]+tmin[i])/86400.
                try:
                    print_fit(fits[part][order][i], args, obsid, bin[i], days, part, order)
                except:
                    pass

def collect_fits(args):

    obsids = args.obsids
    if not obsids:
        obsids = xcal.obsids_years(args.source)[0]

    for i in range(len(obsids)):
        handle_obsid(obsids[i], args)

def main():
    parser = argparse.ArgumentParser(
        description='Collect interleaved calibration observations best fits.'
    )
    parser.add_argument('source', help='mkn421 or pks2155')
    parser.add_argument('det', help='acis or hrc')
    parser.add_argument('-b', '--band', type=int, help='[0-5]')
    parser.add_argument('--wlo', type=float, help='Alternative to --band, short wavelength limit.')
    parser.add_argument('--whi', type=float, help='Alternative to --band, long wavelength limit.')
    parser.add_argument('-o', '--obsids', nargs='+', help='Process only the given obsids')
    args = parser.parse_args()

    if args.band is None and (args.wlo is None or args.whi is None):
        raise ValueError('Must specify either --band or --wlo/--whi')
    if args.band is not None and (args.wlo is not None or args.whi is not None):
        raise ValueError('cannot specify both --band and --wlo/--whi')

    collect_fits(args)

if __name__ == '__main__':
    main()
