import argparse
import numpy as np
import os
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

codes = {
    'ACIS-S' : { 'LEG' : 'l', 'HEG' : 'h',  'MEG' : 'm' },
    'HRC-S' :  { 'LEG' : 's' }
}

# wavelength ranges
# H - 2-18
# H + 2-16 ******
# M - 2-36
# M + 2-31 *****
# L - 2-90
# L + 2-73
# R -
# R +

def qdp_directories():
    return os.environ['QDPPATH'].split(':')

def locate_qdp_directory(obsid):
    qdpdir = None

    for d in qdp_directories():
        qdpneg = '{}/{}_LEG_neg_B0.qdp'.format(d, obsid)
        qdppos = '{}/{}_LEG_pos_B0.qdp'.format(d, obsid)
        if os.path.isfile(qdpneg):
            if qdpdir is not None:
                raise ValueError('Found multiple directories with neg QDP file for obsid {}'.format(obsid))
            else:
                if os.path.isfile(qdppos):
                    qdpdir = d
                else:
                    raise ValueError('Found neg, but not pos, QDP file for obsid {}'.format(obsid))
    if qdpdir is None:
        raise ValueError('Could not find neg QDP file obsid {}'.format(obsid))

    return qdpdir

def qdp_files(obsid, bin):
    qdpdir = locate_qdp_directory(obsid)

    if bin == -1:
        return { o : '{}/{}_LEG_{}_B0.qdp'.format(qdpdir, obsid, o) for o in ('neg', 'pos') }
    else:
        return { o : '{}/{}_LEG_{}_B0_{:02d}.qdp'.format(qdpdir, obsid, o, bin) for o in ('neg', 'pos') }

def qdp_ratios(qdpfile):
    wav, orders = read_qdp(qdpfile)
    ratio = orders[0] / orders.sum(axis=0)
    return wav[::-1], ratio[::-1]

def read_qdp(qdpfile):
    data = np.genfromtxt(qdpfile, skip_header=3, loose=True, unpack=True)
    n = int((data.shape[1]-1)/2)
    return data[0,:n], data[5:,:n]

def get_ratios(obsid, bin):
    files = qdp_files(obsid, bin)
    ratios = { }
    for order in files:
        wav, ratio = qdp_ratios(files[order])
        ratios[order] = { 'wav': wav, 'ratio' : ratio }
    return ratios

def print_header(args):
    print('# source: {}'.format(args.source))
    print('# wavmin: {}'.format(args.wavmin))
    print('# wavmax: {}'.format(args.wavmax))
    print('obsid\tbin\tdays\tinst\torder\tflux\terr')
    print('N\tN\tN\tS\tN\tN\tN')

def handle_obsid(obsid, args):

    global codes

    jnk, pha2 = xcal.pha2_files(obsid) # just binned files
    jnk, evt2 = xcal.evt2_files(obsid) # just binned files
    hdr = util.read_header(jnk)

    detnam, grating = hdr['detnam'], hdr['grating']
    if detnam != 'HRC-S': detnam = 'ACIS-S'

    ratios = None

    arms = ('LEG',)
    if grating == 'HETG':
        arms = ('MEG', 'HEG')

    for arm in arms:

        bin_lo, bin_hi, resp = response.get_response(obsid, arm, maxorder=1)
        bin_lo = bin_lo[0][::-1]
        bin_hi = bin_hi[0][::-1]
        wav = 0.5 * (bin_lo + bin_hi)

        for order in resp: resp[order] = resp[order][0][::-1]

        if args.combine != 1:
            mod = np.mod(bin_lo.size, args.combine)
            if mod:
                bin_lo = bin_lo[:-mod]
                bin_hi = bin_hi[:-mod]
                for order in resp:
                    resp[order] = resp[order][:-mod]

            newshape = (int(bin_lo.size/args.combine), args.combine)

            bin_lo = bin_lo[::args.combine]
            bin_hi = bin_hi[args.combine-1::args.combine]

            for order in resp:
                resp[order] = np.reshape(resp[order], newshape).mean(axis=1)

        orders = { 'neg' : -1, 'pos' : +1 }

        for i in range(len(pha2)):

            if detnam == 'HRC-S':
                ratios = get_ratios(obsid, i)

            t1, t2 = xcal.trange_evt2(evt2[i])
            time = 0.5 * (t1+t2)

            data, hdr = util.read_pha2(pha2[i])

            for order in orders:

                ind = np.where(
                    (data['tg_m']==orders[order]) &
                    (data['tg_part']==util.tg_parts[arm])
                )[0][0]
                src = data['counts'][ind]
                bg = data['background_up'][ind] + data['background_down'][ind]

                src = src[::-1]
                bg = bg[::-1]

                if args.combine != 1:
                    if mod:
                        src = src[:-mod]
                        bg = bg[:-mod]
                    src = np.reshape(src, newshape).sum(axis=1)
                    bg = np.reshape(bg, newshape).sum(axis=1)

                factor = np.ones_like(resp[order])

                if ratios is not None:
                    factor = np.interp(wav, ratios[order]['wav'], ratios[order]['ratio'])
            
                ind = np.where((bin_lo>=args.wavmin) & (bin_hi<args.wavmax))[0]
                rate, rate_err, flux_, flux_err = flux.flux_summed(src[ind], bg[ind], resp[order][ind], bin_lo[ind], bin_hi[ind], hdr, factor[ind])
                print('{}\t{}\t{:.5f}\t{}\t{}\t{:.4g}\t{:.4g}'.format(obsid, i, time/86400., codes[detnam][arm], order, flux_, flux_err))
            #    print('obsid\tbin\tdays\tinst\torder\tflux\terr')



def count_fluxes(args):

    print_header(args)

    obsids, years = xcal.obsids_years(args.source)

    for i in range(len(obsids)):
        handle_obsid(obsids[i], args)

def main():

    parser = argparse.ArgumentParser(
        description='Count up fluxes in interleaved calibration observations.'
    )
    parser.add_argument('-c', '--combine', type=int, default=1, help='Number of spectral bins to group. Default=1')
    parser.add_argument('source', help='mkn421 or pks2155')
    parser.add_argument('wavmin', type=float, help='Minimum wavelength')
    parser.add_argument('wavmax', type=float, help='Maximum wavelength')
    args = parser.parse_args()

    count_fluxes(args)

if __name__ == '__main__':
    main()
