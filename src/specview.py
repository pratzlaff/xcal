import argparse
import numpy as np
import matplotlib.pyplot as plt

import sherpa.utils.err
from sherpa.astro.ui import *
import cxcdm
import ciao_contrib._tools.fileio as fileio

def rows_info(pha2, rows):
    orders = []
    arms = []

    bl = cxcdm.dmTableOpen(pha2)
    nrows = cxcdm.dmTableGetNoRows(bl)
    tg_m = cxcdm.dmTableOpenColumn(bl, 'tg_m')
    tg_part = cxcdm.dmTableOpenColumn(bl, 'tg_part')

    tg_m_vals = cxcdm.dmGetData(tg_m, 1, nrows)
    tg_part_vals = cxcdm.dmGetData(tg_part, 1, nrows)

    orders = [ tg_m_vals[row-1] for row in rows ]
    arms = [ tg_part_vals[row-1] for row in rows ]

    return orders, arms

def main():
    parser = argparse.ArgumentParser(
        description='Plot a spectrum.'
    )
    parser.add_argument('pha2', help='PHA2 file.')
    parser.add_argument('rows', nargs='+', type=int, help='Rows.')
    args = parser.parse_args()

    pha2 = args.pha2
    rows = args.rows

    orders, arms = rows_info(pha2, rows)

    load_pha(pha2, use_errors=True)

    set_analysis('wave')

    for i in range(len(rows)):
        plt.subplot(len(rows), 1, i+1)
        set_default_id(rows[i])
        group_snr(10)
        notice(1.5, 25)
        plot_data(clearwindow=False)

    plt.show()


if __name__ == '__main__':
    main()
