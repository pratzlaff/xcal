import numpy as np
import astropy.io.fits
import argparse
import re
import glob

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from matplotlib import rc, rcParams
rc('text', usetex=True)

datadir='/data/legs/rpete/flight/xcal/data/*'

def evt1a_file(obsid, osip, newgain):
    global datadir

    newgain_str = ''
    if newgain: newgain_str='_newgain'

    osip_str = ''
    if not osip: osip_str='_noosip'

    files = glob.glob(datadir+'/'+obsid+'/tg_reprocess'+newgain_str+'/*_evt1a'+osip_str+'.fits')

    pattern = r'acisf\d+_\d+N\d+_evt1a'+osip_str+'.fits'
    for name in files:
        if re.search(pattern, name): return name
    raise "whoops"

def read_evt1a(evt1a):

    hdulist = astropy.io.fits.open(evt1a)
    data = hdulist['events'].data

    energy = data.field('energy') 
    tg_m = data.field('tg_m')
    tg_mlam = data.field('tg_mlam')
    tg_d = data.field('tg_d')

    hdulist.close()

    return energy, tg_m, tg_mlam, tg_d

def get_orders(obsid, osip, newgain, order_max=None):

    evt1a = evt1a_file(obsid, osip, newgain)
    energy, tg_m, tg_mlam, tg_d = read_evt1a(evt1a)
    order = tg_mlam / (12.39854 / (energy/1e3))

    mask = (energy>0) & (energy==energy) & (tg_mlam==tg_mlam) & (np.abs(tg_d)<0.0008) & (np.abs(order)<5)

    if order_max:
        mask &= np.abs(tg_m)<=order_max

    energy = energy[mask]
    tg_m = tg_m[mask]
    tg_mlam = tg_mlam[mask]
    order = order[mask]

    return order, tg_mlam

def main():

    parser = argparse.ArgumentParser(
        description='FIXME'
    )
    parser.add_argument('-p', '--pdf', help='Save plot to named file.')
    parser.add_argument('-o', '--order_max', type=int, help='Order max.')
    parser.add_argument('obsid')
    args = parser.parse_args()

    orders = {
        'osip' : { },
        'noosip' : { },
        'newgain_osip' : { },
        'newgain_noosip' : { }
    }

    newgain=False
    osip=True
    orders['osip']['order'], orders['osip']['tg_mlam'] = get_orders(args.obsid, osip, newgain, args.order_max)
    osip=False
    orders['noosip']['order'], orders['noosip']['tg_mlam'] = get_orders(args.obsid, osip, newgain, args.order_max)

    newgain=True
    osip=True
    orders['newgain_osip']['order'], orders['newgain_osip']['tg_mlam'] = get_orders(args.obsid, osip, newgain, args.order_max)
    osip=False
    orders['newgain_noosip']['order'], orders['newgain_noosip']['tg_mlam'] = get_orders(args.obsid, osip, newgain, args.order_max)

    ranges = [
        [-48,-45],
        [-43,-40],
        [-27,-24],
        [24,27],
        [45,48]
        ]

    if args.pdf:
        pdf = PdfPages(args.pdf)

    plot_dims = (2, 3)
    fig = plt.figure(figsize = (11, 8.5))

    i = 0
    for r in ranges:
        row = int( i / plot_dims[1] )
        col = i % plot_dims[1]
        plt.subplot2grid(plot_dims, (row, col))

        plt.ylabel(r'$N$')
        plt.xlabel(r'$\textrm{"Order"}$')
        title=r'$\textrm{'+args.obsid+r'}: '+("%.0f" % r[0])+r'<\textrm{TG\_MLAM}<'+("%.0f" % r[1])+'$'
        title=args.obsid+': $'+("%.0f" % r[0])+r' < \textrm{tg\_mlam} < '+("%.0f" % r[1])+'$'
        plt.title(title)

        order, tg_mlam = orders['osip']['order'], orders['osip']['tg_mlam']
        mask = (tg_mlam>r[0]) & (tg_mlam<r[1])
        plt.hist(order[mask], bins=100, histtype='step', label='OSIP', color='r')

        order, tg_mlam = orders['newgain_osip']['order'], orders['newgain_osip']['tg_mlam']
        mask = (tg_mlam>r[0]) & (tg_mlam<r[1])
        plt.hist(order[mask], bins=100, histtype='step', label='OSIP, newgain', color='b')

        # order, tg_mlam = orders['noosip']['order'], orders['noosip']['tg_mlam']
        # mask = (tg_mlam>r[0]) & (tg_mlam<r[1])
        # plt.hist(order[mask], bins=100, histtype='step', label='No OSIP')

        # order, tg_mlam = orders['newgain_noosip']['order'], orders['newgain_noosip']['tg_mlam']
        # mask = (tg_mlam>r[0]) & (tg_mlam<r[1])
        # plt.hist(order[mask], bins=100, histtype='step', label='No OSIP, newgain')

        plt.legend(fontsize=8)
        i += 1

    plt.tight_layout()

    if args.pdf:
        pdf.savefig(fig)
        pdf.close()
    else:
        plt.show()

    plt.close()

    print order.min()
    print order.max()

if __name__ == '__main__':
    main()

