import numpy as np
import argparse
import matplotlib.pyplot as plt

import matplotlib
#matplotlib.rc('text', usetex=True)

import xcal

ANGSTROM, LAMBDA = "Åλ"

def factors_file(args):
    return './results/{}/factors_{}_{}_{}_winc_{:1.0f}.txt'.format(args.caldb, args.type, args.inst, args.source, args.winc)

def main():

    global ANGSTROM

    parser = argparse.ArgumentParser(
        description='Plot factors for the given file'
    )
    parser.add_argument('-o', '--outfile', help='Output file name')
    parser.add_argument('source', help='mkn421|pks2155')
    parser.add_argument('inst', help='Instrument comparison, e.g., sl')
    parser.add_argument('caldb', help='CalDB setup')
    parser.add_argument('type', help='lengths|chi2|subjective|quick|simple|obsflux_quick|obsflux_simple')
    parser.add_argument('--wmin', type=float, help='Columns represent wavlengths ranges, with this as the minimum')
    parser.add_argument('--winc', type=float, help='Columns represent wavlengths ranges, with this as the increment')
    args = parser.parse_args()

    if (args.wmin is not None and args.winc is None) or (args.wmin is None and args.winc is not None):
        raise ValueErro('--wmin/--winc combination is invalid')

    obsids, years = xcal.obsids_years(args.source)
    days = 365.2425*(years-1998)
    groups = xcal.group_days(days)

    # FIXME
    use = {
        'mkn421' : [i for i in range(len(groups))],
        'pks2155' : [i for i in range(len(groups))]
    }
    years = [years[groups[use[args.source][i]][0]] for i in range(len(use[args.source]))]
    print(years)

    m = np.loadtxt(factors_file(args), unpack=True)
    print(m)

    if args.type == 'quick' or args.type == 'simple' or args.type == 'obsflux_quick' or args.type == 'obsflux_simple':
        med = m
        nbands = med.shape[0]
    else:
        min = [m[2*i] for i in range(5)]
        max = [m[2*i+1] for i in range(len(min))]
        med = [0.5*(min[i]+max[i]) for i in range(len(min))]
        nbands = len(med)

    offsets = np.arange(nbands)*0.1 / (nbands-1)
    offsets -= offsets[int(offsets.size/2)]

    if args.winc:
        for i in range(nbands):
            wlo = i * args.winc + args.wmin
            plt.plot(years+offsets[i], med[i], label='{:g}-{:g} {}'.format(wlo, wlo+args.winc, ANGSTROM))
    else:
        for band in 0,2,3,4:
            wlo, whi = xcal.band_wav_limits(band)
            plt.errorbar(years+offsets[band], med[band], 0.5*(max[band]-min[band]), label='{:.1f}-{:.1f} {}'.format(wlo, whi, ANGSTROM))
    plt.xlabel('Year')
    at = { 'lengths' : 'Minimum Length', 'subjective' : 'Subjective Best Match', 'chi2' : r'Minimum $\chi^2$', 'quick' : 'Quick Look Best Match (fit flux)', 'obsflux_quick' : 'Quick Look Best Match (obs flux)', 'simple' : 'Simple Method Best Match (fit flux)', 'obsflux_simple' : 'Simple Method Best Match (obs flux)', }
    plt.ylabel('HRC-S/LETG Scale at {}'.format(at[args.type]))
    plt.title({'pks2155' : 'PKS 2155-304', 'mkn421' : 'Mkn 421'}[args.source])

    plt.legend()
    plt.tight_layout()
    
    if args.outfile:
        plt.savefig(args.outfile)
    else:
        plt.show()
        

if __name__ == '__main__':
    main()
