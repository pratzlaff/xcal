import argparse
import glob
from matplotlib import pyplot as plt
import numpy as np
import os
import re
import sys

import sherpa.utils.err
from sherpa.astro.ui import *
import cxcdm
import ciao_contrib._tools.fileio as fileio

import sherpa_plot
import xcal
import response

# from https://cxc.cfa.harvard.edu/sherpa/threads/grating_hrcsletg/index.html#response

from sherpa.models.model import CompositeModel
from sherpa.astro.instrument import MultiResponseSumModel

def startup_monkey(self, cache):
    pha = self.pha
    if np.iterable(pha.mask):
        pha.notice_response(True)
    self.channel = pha.get_noticed_channels()
    self.mask = pha.get_mask()
    self._get_noticed_energy_list()
    CompositeModel.startup(self, cache)

MultiResponseSumModel.startup = startup_monkey

# Sherpa wstat throws a lot of background warnings
import warnings
warnings.simplefilter('ignore', UserWarning)

TG_PARTS = { 'HEG' : 1, 'MEG' : 2, 'LEG' : 3 }

def nH( src ):
    return {
        'mkn421'  : 0.0161,
        'pks2155' : 0.0142,
        }[src]

def source_def( src, band ):
    models = {
        0 : source_logpar,
    }
    return models.get( band, source_powerlaw )(src)

def source_logpar( src ):
    set_xsabund("wilm")
    set_source(xstbabs.tbabs * xslogpar.lp)
    tbabs.nh = nH( src )
    freeze(tbabs.nh)
    set_par(lp.alpha, val=2, min=0, max=4)
    set_par(lp.beta, val=0, min=-4, max=4)
    set_par(lp.norm, val=0.1, min=0, max=1e24)
    return lp

def source_powerlaw( src ):
    set_xsabund("wilm")
    set_source(xstbabs.tbabs * xspowerlaw.pl)
    tbabs.nh = nH( src )
    freeze(tbabs.nh)
    set_par(pl.PhoIndex, val=2.5, min=-2, max=9)
    set_par(pl.norm, val=0.2, min=0, max=1e24)
    return pl

def read_conf(conffile):
    params = {}
    reading_params = False
    f = open(conffile)
    for line in f:
        if re.search('Param', line):
            reading_params = True
            next(f)
            continue
        if reading_params:
            match = re.search('^\s*(\w+\.\w+)\s+(\S+)\s+(\S+)\s+(\S+)\s*$', line)
            if match:
                param = match.groups()[0]

                try: best = float(match.groups()[1])
                except: best=np.nan

                try: low = float(match.groups()[2])
                except: low=np.nan

                try: high = float(match.groups()[3])
                except: high=np.nan
 
                params[param] = {}
                params[param]['best'] = best
                params[param]['min'] = best + low
                params[param]['max'] = best + high
                
    f.close()
    return params

def get_spec_num(pha2, part, orders):
    global TG_PARTS

    order_int = { 'neg' : -1, 'pos' : 1 }[orders]

    bl = cxcdm.dmTableOpen(pha2)
    nrows = cxcdm.dmTableGetNoRows(bl)    

    spec_num = cxcdm.dmTableOpenColumn(bl, 'spec_num')
    tg_m = cxcdm.dmTableOpenColumn(bl, 'tg_m')
    tg_part = cxcdm.dmTableOpenColumn(bl, 'tg_part')

    spec_num_vals = cxcdm.dmGetData(spec_num, 1, nrows)
    tg_m_vals = cxcdm.dmGetData(tg_m, 1, nrows)
    tg_part_vals = cxcdm.dmGetData(tg_part, 1, nrows)

    return spec_num_vals[(tg_part_vals==TG_PARTS[part]) & (tg_m_vals==order_int)][0]

def pha2_files(dir, args):
    pha2 = glob.glob(dir+'/*_pha2.fits')
    if not args.nobinned:
        pha2.extend(sorted(glob.glob(dir+'/pha2_??.fits')))
    return pha2
#    pha2_binned = sorted(glob.glob(dir+'/pha2_??.fits'))
#    return pha2, pha2_binned

# taken from tgsplit, with modification
def rmf_from_garf(garf):
    grid = fileio.get_keys_from_file(garf)['RESPFILE']
    isgrid = grid.split("(")
    if 2 != len(isgrid) or isgrid[0].lower() != 'grid':
        raise ValueError("Unexpected RESPFILE found in {}".format(infile))
    # strip trailing ")", remove DM filter,
    ff = isgrid[1].split(")")[0].split("[")[0]
    return ff

def locate_garf(obsid, part, tg_m):
    garf = None
    for d in response.arf_directories():
        tmp_garf = glob.glob('{}/{}_{}_{}_garf.fits'.format(d, obsid, part, tg_m))
        if len(tmp_garf):
            if garf is not None:  raise ValueError('Found multiple directories with garf for obsid {}, part {}, tg_m {}'.format(obsid, part, tg_m))
            elif (len(tmp_garf)>1): raise ValueError('Found multiple garfs in directory {}: obsid {}, part {}, tg_m {}'.format(d, obsid, part, tg_m))
            else: garf = tmp_garf[0]
    if garf is None: raise ValueError('Could not find garf for obsid {}, part {}, tg_m {}'.format(obsid, part, tg_m))
    return garf

def garf_files_struct(pha2):
    obsid, obj, detnam, grating, tg_part, tg_m = pha2_info(pha2)

    parts = { 'HETG' : ['MEG', 'HEG'], 'LETG' : ['LEG'] }[grating]

    orders = { 'ACIS' : { 'neg' : [-1], 'pos' : [1] }, 'HRC-' : { 'neg' : tg_m[tg_m<0].tolist(), 'pos' : tg_m[tg_m>0].tolist() } }[detnam[0:4]]
    olist = orders['pos']

    garfs = { }

    for part in parts:
        garfs[part] = { }
        for order in orders:
            garfs[part][order] = { }
            garfs[part][order]['olist'] = olist
            garfs[part][order]['garfs'] = [ locate_garf(obsid, part, tg_m) for tg_m in orders[order] ]
    return garfs

def pha2_obsid(pha2):
    headerkeys = fileio.get_keys_from_file(pha2)
    return headerkeys['OBS_ID']

def pha2_info(pha2):

    headerkeys = fileio.get_keys_from_file(pha2)
    obsid = headerkeys['OBS_ID']
    detnam = headerkeys['DETNAM']
    grating = headerkeys['GRATING']
    obj = headerkeys['OBJECT']

    bl = cxcdm.dmTableOpen(pha2)
    nrows = cxcdm.dmTableGetNoRows(bl)    
    tg_m = cxcdm.dmTableOpenColumn(bl, 'tg_m')
    tg_part = cxcdm.dmTableOpenColumn(bl, 'tg_part')

    tg_m_vals = cxcdm.dmGetData(tg_m, 1, nrows)
    tg_part_vals = cxcdm.dmGetData(tg_part, 1, nrows)

    cxcdm.dmTableClose(bl)

    if re.search(r'^pks\s*2155', obj, re.IGNORECASE):
        obj = 'pks2155'
    elif re.search(r'^mkn\s*421', obj, re.IGNORECASE):
        obj = 'mkn421'
    else:
        raise ValueError("Unrecognized OBJECT={}, from {}".format(obj, pha2))

    return obsid, obj, detnam, grating, tg_part_vals, tg_m_vals

def fitall(args):
    pha2 = pha2_files(args.indir, args)
    obsid = pha2_obsid(pha2[0])

    # FIXME: this stuff doesn't work
    confbases = [ None for p in pha2 ]
    try:
        base, base_bins = xcal.fit_basenames(obsid, args.part, args.orders, args.band, args.wlo, args.whi)
        bases = [ base ]
        if not args.nobinned:
            bases.extend(base_bins)
        if len(bases) == len(pha2): confbases=bases
    except:
        pass

    for i in range(len(pha2)):
        fitone(pha2[i], args, confbase=confbases[i])

def fitone(pha2, args, confbase=None):

    obsid, obj, detnam, grating, tg_part_vals, tg_m_vals = pha2_info(pha2)
    garfs_struct = garf_files_struct(pha2)

    binspec=''
    match = re.match( r'.*/pha2_(\d+).fits', pha2 )
    if match:
        binspec = '_{}'.format(match.groups()[0])

    for part in (args.part,):
    #for part in garfs_struct:

        for orders in (args.orders,):
        #for orders in garfs_struct[part]:

            if args.band is not None:
                outbase = '{}_{}_{}_B{}{}'.format(obsid, part, orders, args.band, binspec)
            else:
                outbase = '{}_{}_{}_w{:g}-{:g}{}'.format(obsid, part, orders, args.wlo, args.whi, binspec)
                

            pdffile = args.pdir + '/' + outbase + '.pdf'
            conffile = args.fdir + '/' + outbase + '.conf'
            paramsfile = args.fdir + '/' + outbase + '.params'
            fluxfile = args.fdir + '/' + outbase + '.flux'
            fluxesfile = args.fdir + '/' + outbase + '_fsamp.npy'
            show_allfile = args.fdir + '/' + outbase + '.show_all'
            save_allfile = args.fdir + '/' + outbase + '.save_all'

            garfs = garfs_struct[part][orders]['garfs']
            olist = garfs_struct[part][orders]['olist']
            rmfs = [ rmf_from_garf(arf) for arf in garfs]

            load_pha(1, pha2, use_errors=True)
            set_default_id(get_spec_num(pha2, part, orders))

            group_counts(1)

            set_stat('wstat')

            load_multi_arfs(garfs, olist)
            load_multi_rmfs(rmfs, olist)

            set_analysis('wave')
            ignore()
            if args.band is not None:
                wlo, whi = xcal.band_wav_limits(args.band)
            else:
                wlo, whi = args.wlo, args.whi
            notice(wlo, whi)

            # unabsorbed component
            unabs_comp = source_def(obj, args.band)

            if (args.usebase):
                try:
                    basename = xcal.fit_basename(obsid, args.part, args.orders, args.band, args.wlo, args.whi)
                    initfile = basename + '.conf'
                    initdict = read_conf(initfile)
                    factor=1.0
                    print(initdict)
                    for param in initdict:
                        val = initdict[param]['best']
                        min = initdict[param]['best'] - factor*(initdict[param]['best']-initdict[param]['min'])
                        max = initdict[param]['best'] + factor*(initdict[param]['max']-initdict[param]['best'])
                        print("SETTING {}: val={}, min={}, max={}".format(param, val, min, max))
                        set_par(eval(param), val=val, min=min, max=max)
                except ValueError:
                    pass

            # if confbase is not None:
            #     try:
            #         confdict = xcal.read_fit(confbase)
            #         print("CONFDICT=", condict)
            #         # FIXME
            #         params = { 'alpha' : lp.alpha, 'beta' : lp.beta, 'norm' : lp.norm }
            #         for param in confdict:
            #             if param == 'flux': continue
            #             set_par(params[param], val=confdict[param]['best'])
            #     except:
            #         pass

            if not args.nofit:
                for f in paramsfile, fluxfile, show_allfile, save_allfile:
                    try:
                        os.remove(f)
                    except FileNotFoundError:
                        pass

                fit()
                show_fit(outfile=paramsfile, clobber=True)

                flux = calc_energy_flux(wlo, whi)
                f=open(fluxfile, 'w' )
                f.write("energy flux = %s\n" % flux)
                f.close()

                if not args.noconf:
                    try:
                        os.remove(conffile)
                    except FileNotFoundError:
                        pass
                    set_conf_opt("max_rstat", 4)
                    try:
                        conf()
                        show_conf( outfile=conffile, clobber=True )
                    except sherpa.utils.err.EstErr:
                        pass

                if not args.nosamp:
                    try:
                        os.remove(fluxesfile)
                    except FileNotFoundError:
                        pass
                    s = sample_flux(unabs_comp, wlo, whi, num=1000)
                    fluxes = np.sort(s[2][:,0])
                    np.save(fluxesfile, fluxes)

                show_all( outfile=show_allfile, clobber=True )
                #save_all( outfile=save_allfile, clobber=True )

            if not args.noplot:
                try:
                    os.remove(pdffile)
                except FileNotFoundError:
                    pass
                group_snr(args.snr)
                title = obsid+': '+{'ACIS':'ACIS-S','HRC-':'HRC-S'}[detnam[0:4]]+'/'+args.part+', TG_M='+{'neg':'-', 'pos':'+'}[args.orders]+'1'
                sherpa_plot.plot_range(wlo, whi, title, units='wavelength')
                plt.savefig(pdffile)
                ungroup()

def main():
    
    parser = argparse.ArgumentParser(
        description='Perform Sherpa fits, given an input directory.'
    )
    parser.add_argument('-p', '--pdir', default='./plots', help='Directory to which plots are saved')
    parser.add_argument('-f', '--fdir', default='./fits', help='Directory to which fit information are saved')
    parser.add_argument('--nofit', action='store_true', help='Do not perform a fit.')
    parser.add_argument('--nobinned', action='store_true', help='Do not perform binned fits.')
    parser.add_argument('--noconf', action='store_true', help='Do not calculate confidence interval.')
    parser.add_argument('--nosamp', action='store_true', help='Do not sample fluxes.')
    parser.add_argument('-u', '--usebase', action='store_true', help='For binned fits, set initial parameters from the base best fit.')
    parser.add_argument('--noplot', action='store_true', help='Do not write a plot.')
    parser.add_argument('--snr', default=20., type =float, help='SNR grouping value for plot. Default is 5.')
    parser.add_argument('indir', help='Directory containing reprocessed data.')
    parser.add_argument('part', help='HEG, MEG, or LEG.')
    parser.add_argument('orders', help='neg or pos.')
    parser.add_argument('-b', '--band', type=int, help='[0-5]')
    parser.add_argument('--wlo', type=float, help='Alternative to --band, short wavelength limit.')
    parser.add_argument('--whi', type=float, help='Alternative to --band, long wavelength limit.')
    args = parser.parse_args()

    if args.band is None and (args.wlo is None or args.whi is None):
        raise ValueError('Must specify either --band or --wlo/--whi')
    if args.band is not None and (args.wlo is not None or args.whi is not None):
        raise ValueError('cannot specify both --band and --wlo/--whi')

    xcal.mkdir_p(args.fdir)
    xcal.mkdir_p(args.pdir)

    os.environ['FITPATH'] = args.fdir # used to look for previous fit

    fitall(args)

    # pha2, pha2_binned = pha2_files(args.indir)

    # garfs_struct = garf_files_struct(pha2)
    # print garfs_struct

    # for part in garfs_struct:
    #     for orders in garfs_struct[part]:
    #         arfs = garfs_struct[part][orders]
    #         rmfs = [ rmf_from_garf(arf) for arf in arfs]
    #         print arfs
    #         print rmfs
    #         print get_spec_num(pha2, part, orders)

if __name__ == '__main__':
    main()
