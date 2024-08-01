import argparse
import numpy as np
import matplotlib.pyplot as plt

import response

def plot_areas():
    
    parser = argparse.ArgumentParser(
        description='Plot areas for 20 Year poster.'
    )
    parser.add_argument('-o', '--outfile', help='Output file name.')
    args = parser.parse_args()

    # PKS 2155-304, ACIS-S/LETG, 2000-05-31
    aw1_2000, aw2_2000, aa_2000 = response.read_arf('./arfs/acis/1703_LEG_1_garf.fits')
    ae_2000 = 12.39854/0.5/(aw1_2000+aw2_2000)

    # PKS 2155-304, HRC-S/LETG, 2000-05-31
    hw1_2000, hw2_2000, ha_2000 = response.read_arf('./arfs/qe_N0014_qeu_N0013/1704_LEG_1_garf.fits')
    he_2000 = 12.39854/0.5/(hw1_2000+hw2_2000)

    # Mkn 421, ACIS-S/LETG, 2019-07-05
    aw1_2019, aw2_2019, aa_2019 = response.read_arf('./arfs/acis/21817_LEG_1_garf.fits')
    ae_2019 = 12.39854/0.5/(aw1_2019+aw2_2019)

    # Mkn 421, HRC-S/LETG, 2019-07-05
    hw1_2019, hw2_2019, ha_2019 = response.read_arf('./arfs/qe_N0014_qeu_N0013/21818_LEG_1_garf.fits')
    he_2019 = 12.39854/0.5/(hw1_2019+hw2_2019)

    #plt.title('Effective Areas, LETG Positive 1st Order')
    plt.xlabel('Energy (keV)')
    plt.ylabel(r'Area ($\mathregular{cm^2}$)')

    plt.xscale('log')
    plt.yscale('log')
    plt.ylim(1, 100)

    plt.plot(ae_2000, aa_2000, 'k--', label='ACIS-S, 2000')
    plt.plot(ae_2019, aa_2019, 'k-', label='ACIS-S, 2019')

    plt.plot(he_2000, ha_2000, 'r--', label='HRC-S, 2000')
    plt.plot(he_2019, ha_2019, 'r-', label='HRC-S, 2019')

    plt.legend()

    plt.tight_layout()

    if args.outfile:
        plt.savefig(args.outfile)
    else:
        plt.show()

if __name__ == '__main__':
    plot_areas()
