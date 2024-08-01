import argparse
import numpy as np
import matplotlib.pyplot as plt

import util
import xcal

def main():

    parser = argparse.ArgumentParser(
        description='Plot plus/minus count ratio vs wavelength for each epoch.'
    )
    parser.add_argument('-o', '--outfile', help='Output file name.')
    parser.add_argument('-c', '--combine', type=int, default=128, help='Bins to combine.')
    parser.add_argument('--wmin', default=2., help='Minimum wavelength.')
    parser.add_argument('--wmax', default=24., help='Maximum wavelength.')
    parser.add_argument('src', help='mkn421|pks2155')
    args = parser.parse_args()
    

    # get obsids

    obsids, years = xcal.obsids_years(args.src)
    groups = xcal.group_days((years-1998)*365.2425)

    for i in range(len(groups)):
        date_obs = None
        group = groups[i]
        wav, counts = None, {}
        for j in range(len(group)):
            obsid = obsids[group[j]]
            pha2 = xcal.pha2_files(obsid)[0]
            detnam, grating = xcal.instruments_pha2(pha2)
            if detnam != 'HRC-S':
                continue
            data, hdr = util.read_pha2(pha2)
            if date_obs is None:
                date_obs = hdr['date-obs'][:10]
            orders = {
                'pos' : 1,
                'neg' : -1,
            }
            for order in orders:
                ind = np.where((data['tg_m']==orders[order]) &
                               (data['tg_part']==3)
                )[0][0]
                src = data['counts'][ind]
                bg = (data['background_up'][ind]+data['background_down'][ind])/(hdr['backscup']+hdr['backscdn'])
                if wav is None:
                    wav = 0.5 * (data['bin_lo'][ind] + data['bin_hi'][ind])
                if order not in counts:
                    counts[order] = {
                        'src' : src,
                        'bg' : bg
                    }
                else:
                    counts[order]['src'] += src
                    counts[order]['bg'] += bg

        mask = (wav>args.wmin) & (wav<args.wmax)
        wav = wav[mask]
        for order in orders:
            for type in counts[order]:
                counts[order][type] = counts[order][type][mask]

        if args.combine != 1:
            newshape = (int(wav.size / args.combine), args.combine)
            mod = np.mod(wav.size, args.combine)
            if mod:
                wav = wav[:-mod]
            wav = np.reshape(wav, newshape).mean(axis=1)

            for order in orders:
                if mod:
                    for type in counts[order]:
                        counts[order][type] = counts[order][type][:-mod]

                counts[order]['src'] = np.reshape(counts[order]['src'], newshape).sum(axis=1)
                counts[order]['bg'] = np.reshape(counts[order]['bg'], newshape).sum(axis=1)

        net, net_err = {}, {}
        for order in orders:
            net[order] = counts[order]['src'] - counts[order]['bg']
            net_err[order] = np.sqrt(counts[order]['src'] + counts[order]['bg'])
        ratio = net['pos'] / net['neg']
        ratio_err = ratio * np.sqrt( (net_err['pos']/net['pos'])**2 + (net_err['neg']/net['neg'])**2 )
        plt.errorbar(wav, ratio, ratio_err, label=date_obs)

    sources = {
        'mkn421' : 'Mkn 421',
        'pks2155' : 'PKS 2155-304'
    }

    ANGSTROM, LAMBDA = "Åλ"

    plt.title('{}: HRC-S/LETG +/- Orders Count Ratios'.format(sources[args.src]))
    plt.xlabel('{} ({})'.format(LAMBDA, ANGSTROM))
    plt.ylabel('Count Ratio: Plus / Minus Orders')
    plt.legend()
    plt.tight_layout()

    if args.outfile:
        plt.savefig(args.outfile)
    else:
        plt.show()

if __name__ == '__main__':
    main()
