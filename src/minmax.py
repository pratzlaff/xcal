import argparse
import numpy as np

def minmax_file(args):
    return './minmax/{}_{}_{}'.format(args.src,
                                      args.type,
                                      args.caldb,
    )

def minmax(args):
    d = np.loadtxt(minmax_file(args), unpack=True)
    mincol = 2*(args.band)
    return d[mincol][args.epoch-1], d[mincol+1][args.epoch-1]

def main():
    parser = argparse.ArgumentParser(
        description='Print out scale length/plot minmax.'
    )
    parser.add_argument('src', help='mkn421|pks2155')
    parser.add_argument('type', help='lengths|plots')
    parser.add_argument('caldb', help='qe_N0014_qeu_N0013')
    parser.add_argument('band', type=int, help='[0-5]')
    parser.add_argument('epoch', type=int, help='[1-9]')
    args = parser.parse_args()

    min, max = minmax(args)
    print('{:g} {:g}'.format(min, max))

if __name__ == '__main__':
    main()
