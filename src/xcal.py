import re
import os, errno, sys
import numpy as np
import glob

try:
    import astropy.io.fits
except:
    pass

import util

gsep = 75/86.4 # minimum number of days between separate groups
basedir = '/data/legs/rpete/flight/xcal'
datadir = basedir + '/data'

util.datadir = datadir

configs = {
    'ACIS' : { 'LEG' : 'l',
                 'MEG' : 'm',
                 'HEG' : 'h'
                 },
    'HRC' : { 'LEG' : 's' },
    's' : { 'HRC' : 'LEG' },
    'l' : { 'ACIS' : 'LEG' },
    'm' : { 'ACIS' : 'MEG' },
    'h' : { 'ACIS' : 'HEG' }
}

def e2w(energy):
    return 12.39854 / energy

def w2e(wav):
    return e2w(wav)

def band_energy_limits(band):
    return {
        0 : (.5, 8.0),    # 1.5-25 AA
        1 : (0.33, 0.54), #  23-38 AA, a bit of MEG, but zero HEG
        2 : (0.54, 0.85), #  15-23 AA, still no substantive HEG
        3 : (0.85, 1.5),  #   8-15 AA, full HEG coverage
        4 : (1.5, 4.0),   #    3-8 AA
        5 : (4.0, 10.0),  #  1.2-3 AA, not many counts anywhere
        # 6 : (.26, 0.514)
        }[band]

def band_wav_limits(band):
    return [ e2w(e) for e in band_energy_limits(band) ][::-1]

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise

def fit_directories():
    return os.environ['FITPATH'].split(':')

def locate_fit_directory(obsid, part, orders, band, wlo, whi):
    confdir = None

    for d in fit_directories():
        if band is not None:
            globstr = '{}/{}_{}_{}_B{}.conf'.format(d, obsid, part, orders, band)
        else:
            globstr = '{}/{}_{}_{}_w{:g}-{:g}.conf'.format(d, obsid, part, orders, wlo, whi)
        conf = glob.glob(globstr)
        if len(conf):
            if confdir is not None:  raise ValueError('Found multiple directories with fit conf for obsid {}'.format(obsid))
            elif (len(conf)>1): raise ValueError('Found multiple fit confs in directory {}: obsid {}'.format(d, obsid))
            else: confdir = d
    if confdir is None: raise ValueError('Could not find fit conf for obsid {}, part {}, orders {}, band {}, wlo {:g}, whi {:g} - last globstr is {}'.format(obsid, part, orders, band, wlo, whi, globstr))
    return confdir

def generic_files(obsid, type):
    global datadir
    pha2 = glob.glob('{}/*/{}/tg_reprocess/*_{}.fits'.format(datadir, obsid, type))[0]
    binned = sorted(glob.glob('{}/*/{}/tg_reprocess/{}_??.fits'.format(datadir, obsid, type)))
    return pha2, binned

def pha2_files(obsid):
    return generic_files(obsid, 'pha2')

def evt2_files(obsid):
    return generic_files(obsid, 'evt2')

def get_sources():
    return ('mkn421', 'pks2155')

def ngroups(src):
    obsids, years = obsids_years(src)
    return len(group_days(365.2425*(years-1998)))

def obsids_years(source):

    if 'data' not in obsids_years.__dict__:
        obsids_years.data = { }

    if source not in obsids_years.data:

        obsids, years = [], []

        for obsid in get_obsids(source):
            pha2 = util.pha2_file(obsid)
            hdr = util.read_header(pha2)
            obsids.append(obsid)
            years.append(hdr['year'])

        obsids = np.array(obsids)
        years = np.array(years)

        si = np.argsort(years)

        obsids_years.data[source] = { 'obsids' : obsids[si], 'years' : years[si] }


    return obsids_years.data[source]['obsids'], obsids_years.data[source]['years']

def get_obsids(source):
    global basedir

    ofile = '{}/obsids/{}'.format(basedir, source)
    obsids = np.loadtxt(ofile, usecols=(0,), dtype=int)
    return obsids#[57:61]

def trange_evt2(evt2):
    hdulist = astropy.io.fits.open(evt2)
    time = hdulist['events'].data.field('time')
    hdulist.close()
    return time[0], time[-1]

def trange_obsid(obsid):
    evt2_main, evt2_binned = evt2_files(obsid)
    tmin, tmax = [], []
    for f in evt2_binned:
        tmin_, tmax_ = trange_evt2(f)
        tmin.append(tmin_)
        tmax.append(tmax_)
    return tmin[0], tmax[-1], tmin, tmax

def obsids_tranges(source):
    obsids = get_obsids(source)
    tmin, tmax, tmin_bins, tmax_bins = [], [], [], []

    for obsid in obsids:
        tmin_, tmax_, tmin_bins_, tmax_bins_ = trange_obsid(obsid)
        tmin.append(tmin_)
        tmax.append(tmax_)
        tmax_bins.append(tmin_bins_)
        tmin_bins.append(tmax_bins_)

    sorti = sorted(range(len(tmin)), key=lambda ix: tmin[ix])

    obsids = [obsids[i] for i in sorti]
    tmin = [tmin[i] for i in sorti]
    tmax = [tmax[i] for i in sorti]
    tmin_bins = [tmin_bins[i] for i in sorti]
    tmax_bins = [tmax_bins[i] for i in sorti]

    return obsids, tmin, tmax, tmin_bins, tmax_bins

def date_obs_obsid(obsid):
    pha2, jnk = pha2_files(obsid)
    hdulist = astropy.io.fits.open(pha2)
    date_obs = hdulist['spectrum'].header['DATE-OBS']
    hdulist.close()
    return date_obs

def instruments_pha2(pha2):
    hdulist = astropy.io.fits.open(pha2)
    header = hdulist['spectrum'].header
    detnam = header['detnam']
    grating = header['grating']
    hdulist.close()

    return detnam, grating

def instruments_obsid(obsid):
    pha2, jnk = pha2_files(obsid)
    return instruments_pha2(pha2)

def interleaved_obsids(obsids):
    detector = None
    for obsid in obsids:
        detnam, grating = instruments_obsid(obsid)
        if detector is None: detector = detnam[0:4]
        else:
            if detnam[0:4] != detector: return True
    return False

def group_days(days):
    global gsep
    sep = 3
    gaps = np.where((days[1:]-days[:-1])>gsep)[0]
    splits = np.split(np.arange(days.size), gaps+1)
    return splits

    ii = np.where((days[1:]-days[:-1])>gsep)[0]

    groups = [np.arange(ii[0]+1)]

    for i in range(1,ii.size-1):
        groups.append(np.arange(ii[i]-ii[i-1])+ii[i-1]+1)
    groups.append(np.arange(days.size-1-ii[-1])+ii[-1]+1)
    return groups

def group_obsids_tmins(obsids, tmins):
    global gsep
    groups = []

    this_group = [0]
    for i in range(1,len(obsids)):

        # if more than three days past the last observation, start over
        if (tmins[i]-tmins[this_group[-1]]) > gsep*86400:
            groups.append(this_group)
            this_group = [i]
        else:
            this_group.append(i)

    groups.append(this_group)

    return groups

def read_flux(fluxfile):
    f = open(fluxfile)
    line = f.readline()
    f.close()
    flux = float(re.search('energy flux = (.*)', line).groups()[0])
    return flux

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
            match = re.search('^\s*\w+\.(\w+)\s+(\S+)\s+(\S+)\s+(\S+)\s*$', line)
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

def fit_basename(obsid, part, orders, band, wlo, whi):
    fdir = locate_fit_directory(obsid, part, orders, band, wlo, whi)
    if band is not None:
        return '{}/{}_{}_{}_B{}'.format(fdir, obsid, part, orders, band)
    else:
        return '{}/{}_{}_{}_w{:g}-{:g}'.format(fdir, obsid, part, orders, wlo, whi)

def fit_basenames(obsid, part, orders, band, wlo, whi):
    basename = fit_basename(obsid, part, orders, band, wlo, whi)
    basename_bins = [ b[:-5] for b in sorted(glob.glob(basename + '_??.conf')) ]
    return basename, basename_bins

def read_fsamp(fsamp):
    # FIXME: could maybe not assume flux samples are sorted, and be a
    # bit more rigorous about the median
    fluxes = np.load(fsamp)
    offset=int(fluxes.size*(1-.6827)/2)
    return  fluxes[int(fluxes.size/2)], fluxes[offset], fluxes[-offset]

def read_fit(basename):
    conffile = basename + '.conf'
    fluxfile = basename + '.flux'
    fsampfile = basename + '_fsamp.npy'

    conf = read_conf(conffile)

    try:
        med_, min_,  max_ = read_fsamp(fsampfile)
        conf['flux'] = { 'best' : med_, 'min' : min_, 'max' : max_ }
    except:
        sys.stderr.write("Could not read flux samples file '{}'\n".format(fsampfile))
        flux = read_flux(fluxfile)
        conf['flux'] = { 'best' : flux, 'min' : np.nan, 'max' : np.nan }
        try: conf['flux']['min'] = flux * conf['norm']['min'] / conf['norm']['best']
        except: pass
        try: conf['flux']['max'] = flux * conf['norm']['max'] / conf['norm']['best']
        except: pass

    return conf

def get_fits(obsid, part, orders, band, wlo, whi):

    basename, basename_bins = fit_basenames(obsid, part, orders, band, wlo, whi)

    fit_main = read_fit(basename)
    fit_bins = [read_fit(b) for b in basename_bins]
    return fit_main, fit_bins

def get_fits_obsids(obsids, tmins, tmaxs, band, wlo, whi):

    """
return a structure of the form
    { 'ACIS' :
      {
          'LEG' :
          {
              'tmin' : [], 'tmax' : [],
              'neg' :
              {
                  'alpha' : { 'best' : [], 'min' : [], 'max' : []},
                  'beta' : { 'best' : [], 'min' : [], 'max' : []},
                  'norm' : { 'best' : [], 'min' : [], 'max' : []}
              }
              'pos' :
              {
                  'alpha' : { 'best' : [], 'min' : [], 'max' : []},
                  'beta' : { 'best' : [], 'min' : [], 'max' : []},
                  'norm' : { 'best' : [], 'min' : [], 'max' : []}
              }
          }
          'MEG' : { ... }
          'HEG' : { ... }
      }
      'HRC' : { 'LEG' : { ... } }
    }

    """

    data = { }

    for i in range(len(obsids)):

        obsid = obsids[i]
        tmin = tmins[i]
        tmax = tmaxs[i]
        
        detnam, grating = instruments_obsid(obsid)
        det = { 'ACIS' : 'ACIS', 'HRC-' : 'HRC' }[detnam[0:4]]
        parts = { 'LETG' : ['LEG'], 'HETG' : ['MEG', 'HEG'] }[grating]

        for part in parts:

            tlimits_extended = False

            for orders in 'neg', 'pos':
                try:
                    fit_main, fit_bins = get_fits(obsid, part, orders, band, wlo, whi)
                except:
                    continue

                if det not in data: data[det] = { }
                if part not in data[det]:
                    data[det][part] = {
                        'tmin' : [], 'tmax' : [],
                        'neg' : {}, 'pos' : {}
                    }

                if not tlimits_extended:
                    data[det][part]['tmin'].extend(tmin)
                    data[det][part]['tmax'].extend(tmax)
                    tlimits_extended = True

                for fit in fit_bins:

                    for param in fit: # e.g., 'alpha', 'beta', 'norm', 'flux'

                        if param not in data[det][part][orders]:
                            data[det][part][orders][param] = { }

                        for key in fit[param]: # 'best', 'min', 'max'

                            if key not in data[det][part][orders][param]:
                                data[det][part][orders][param][key] = []

                            data[det][part][orders][param][key].append(fit[param][key])
        

    return data

def collect_fits(source, band, wlo, whi):

    # lists of obsids and time ranges
    obsids, tmin, tmax, tmin_bins, tmax_bins = obsids_tranges(source)

    # split into interleaved groups
    groups = group_obsids_tmins(obsids, tmin)

    o_ret = [] # list of lists of obsids
    f_ret = [] # list of lists of fits

    for g in groups:
        obsids_ = [obsids[i] for i in g]

        if not interleaved_obsids(obsids_): continue

        tmin_ = [tmin[i] for i in g]
        tmax_ = [tmax[i] for i in g]
        tmin_bins_ = [tmin_bins[i] for i in g]
        tmax_bins_ = [tmax_bins[i] for i in g]

        fits = get_fits_obsids(obsids_, tmin_bins_, tmax_bins_, band, wlo, whi)

        o_ret.append(obsids_)
        f_ret.append(fits)

    return o_ret, f_ret

def rdb_hdr_setwav(hdr):
    if 'wlo' in hdr:
        return
    if 'band' in hdr:
        wlo, whi = xcal.band_wav_limits(int(hdr['band']))
    elif 'wavmin' in hdr:
        wlo, whi = hdr['wavmin'], hdr['wavmax']
    else:
        raise ValueError('no appropriate keyword found in header: ' + hdr)

    hdr['wlo'] = wlo
    hdr['whi'] = whi

def read_obsflux_rdb(rdbfile):
    hdr = {}
    fh = open(rdbfile)
    try:
        while True:
            line = fh.readline()
            if not line: break
            line = line.rstrip()
            if re.match(r'^#', line):
                match = re.match(r'^#\s+(\w+):\s+(\S+)', line)
                if match:
                    hdr[match.groups()[0]] = match.groups()[1]
            else:
                cols = line.split('\t')
                types = fh.readline().strip().split('\t')
                dtype='int,int,f8,U1,U3,f8,f8'
                #data = np.genfromtxt(fh, dtype=dtype, names=cols, unpack=True)
                data = np.loadtxt(fh, dtype=dtype, unpack=True)
    except:
        raise
    finally:
        fh.close()

    rdb_hdr_setwav(hdr)

    data = { cols[i] : data[i] for i in range(len(cols)) }

    # sequence time
    sorti = np.argsort(data['days'])
    for k in data:
        data[k] = data[k][sorti]

    # FIXME: make retention of full fit an option?
    ii = (data['bin'] >= 0)
    for k in data:
        data[k] = data[k][ii]

    return hdr, data

def read_fits_rdb(rdbfile, dtype=None):
    hdr = {}
    fh = open(rdbfile)
    try:
        while True:
            line = fh.readline()
            if not line: break
            line = line.rstrip()
            if re.match(r'^#', line):
                match = re.match(r'^#\s+(\w+):\s+(\S+)', line)
                if match:
                    hdr[match.groups()[0]] = match.groups()[1]
            else:
                cols = line.split('\t')
                types = fh.readline().strip().split('\t')
                hdr['params'] = cols[5::3]
                if dtype is None:
                    dtype='int,int,f8,U1,U3,' + ','.join(['f8']*(len(cols)-5))
                #data = np.genfromtxt(fh, dtype=dtype, names=cols, unpack=True)
                data = np.loadtxt(fh, dtype=dtype, unpack=True)
    except:
        return read_fits_rdb(rdbfile, dtype='int,int,f8,U1,U3,f8,f8')

    finally:
        fh.close()

    rdb_hdr_setwav(hdr)

    data = { cols[i] : data[i] for i in range(len(cols)) }

    # sequence time
    sorti = np.argsort(data['days'])
    for k in data:
        data[k] = data[k][sorti]

    # FIXME: make retention of full fit an option?
    ii = (data['bin'] >= 0)
    for k in data:
        data[k] = data[k][ii]

    return hdr, data

def read_fits_rdbs(rdbs):
    hdr, data = read_fits_rdb(rdbs[0])

    for i in range(1, len(rdbs)):
        jnk, d = read_fits_rdb(rdbs[i])
        for k in data:
            data[k] = np.append(data[k],d[k])

    # sequence time
    ii = np.argsort(data['days'])
    for k in data:
        data[k] = data[k][ii]

    return hdr, data
