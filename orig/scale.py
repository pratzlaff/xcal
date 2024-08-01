import numpy as np
import os, errno
import glob
import astropy.io.fits

datadir='/data/legs/rpete/flight/xcal/data/*'

HEG=1
MEG=2
LEG=3

def band_range(band):
    return {
        0: (0.5, 8),
        1: (0.33, 0.54),
        2: (0.54, 85),
        3: (0.85, 1.5),
        4: (1.5, 4.0),
        5: (4.0, 10.0)
        }[band]

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise

def get_gaps(time):

    gaps, = np.where(time[1:] - time[0:-1] > 10)
    splits = np.split(time, gaps+1)

    start = np.zeros(gaps.size+1, dtype=np.long)
    stop = start.copy()

    start[1:] = gaps+1
    stop[0:-1] = gaps ; stop[-1] = time.size-1

    return start, stop

def readflux(file):
    obsid, days, pos, pos_low, pos_high, neg, neg_low, neg_high = np.loadtxt(file,
        skiprows=10, unpack=True )
    flux = pos + neg
    flux_err = 0.5 * np.sqrt( (pos_high-pos_low)**2 + (neg_high-neg_low)**2 )

    return { 'obsid':obsid, 'time':days, 'flux':flux, 'err':flux_err }

def get_sim_times(data, n):

    time = np.concatenate([data[det]['time'] for det in data])
    sorti = np.argsort(time)
    time = time[sorti]

    start, stop = get_gaps(time)

    sim_times = np.zeros((start.size, n))

    for i in range(start.size):
        time_i = time[start[i]:stop[i]+1]
        tmin = time_i.min()
        tmax = time_i.max()
        span = tmax-tmin
        sim_times[i]= np.arange(n) * span / n + tmin

    return sim_times

def pha2_file(obsid):
    global datadir
    return glob.glob(datadir+'/'+obsid+'/tg_reprocess/*_pha2.fits')[0]

def resp_files(obsid, grating):
    global datadir
    global HEG, MEG, LEG

    grating_string = { HEG : 'HEG', MEG : 'MEG', LEG : 'LEG' }[grating]
    return {
        'neg' : glob.glob(datadir+'/'+obsid+'/tg_reprocess/*'+grating_string+'_-1_garf.fits')[0],
        'pos' : glob.glob(datadir+'/'+obsid+'/tg_reprocess/*'+grating_string+'_1_garf.fits')[0]
        }

def get_resp(obsid, grating):

    files = resp_files(obsid, grating)

    retval = { 'neg':{}, 'pos':{} }
    for order in files:
        hdulist = astropy.io.fits.open(files[order])
        header = hdulist['specresp'].header
        data = hdulist['specresp'].data

        for field in 'specresp', 'bin_lo', 'bin_hi':
            retval[order][field] = data.field(field)

        hdulist.close()

    return retval

def get_spectra(obsids, grating, combine):

    global HEG, MEG, LEG

    spectra = {'neg':{}, 'pos':{}}
    wav = None

    orders = {'neg':-1, 'pos':+1}

    for i in xrange(len(obsids)):

        obsid = str(obsids[i])

        pha2 = pha2_file(obsid)

        resp = get_resp(obsid, grating)

        hdulist = astropy.io.fits.open(pha2)
        header = hdulist['spectrum'].header
        data = hdulist['spectrum'].data

        rows = {}
        for order in spectra:
            rows[order] = np.where(
                np.logical_and(
                    data.field('tg_m')==orders[order],
                    data.field('tg_part')==grating)
            )[0][0]

        bg = header['backscal'] / (header['backscup']+header['backscdn']) * \
             (data.field('background_up') + data.field('background_down') )

        s = { 'neg':{}, 'pos':{} }

        for order in spectra:
            s[order]['bg'] = bg[rows[order]]
            for field in 'counts', 'bin_lo', 'bin_hi':
                s[order][field] = data.field(field)[rows[order]]

            s[order]['rate'] = (s[order]['counts'] - s[order]['bg']) \
                               / header['livetime']

            # get arf, interpolate onto spectra grid, calculate fluxes

            energy = np.ma.masked_invalid(
                bins_to_energy(s[order]['bin_lo'],
                               s[order]['bin_hi']) )

            resp_energy = bins_to_energy(resp[order]['bin_lo'],
                                         resp[order]['bin_hi'])

            s[order]['specresp'] = np.interp(energy,
                                             resp_energy,
                                             resp[order]['specresp'])

            s[order]['flux'] = s[order]['rate'] * energy * 1.602e-9 / \
                               s[order]['specresp']

            s[order]['lambda'] = 0.5 * (s[order]['bin_hi']+s[order]['bin_lo'])

            # rebin everything
            newshape = (s[order]['flux'].shape[0]/combine, combine)

            for key in 'flux', 'lambda', 'rate', 'counts', 'bg':
                s[order][key] = np.reshape(s[order][key], newshape).sum(axis=1)

            s[order]['lambda'] /= newshape[1]
            s[order]['rate_err'] = np.sqrt(s[order]['counts'] + s[order]['bg']) \
                                   / header['livetime']
            s[order]['flux_err'] = s[order]['flux'] * s[order]['rate_err'] \
                                   / s[order]['rate']

            delta_lambda = np.abs(
                np.median(s[order]['lambda'][:-1] - s[order]['lambda'][1:]) )

            for key in 'flux', 'flux_err':
                s[order][key] /= delta_lambda

            if not spectra[order].keys():
                wav = s[order]['lambda']
                zeros = np.zeros((len(obsids), s[order]['flux'].size))
                for key in 'rate', 'rate_err', 'flux', 'flux_err':
                    spectra[order][key] = np.copy(zeros)
                spectra[order]['livetime'] = np.zeros(len(obsids))


            spectra[order]['livetime'][i] = header['livetime']
            for key in 'rate', 'rate_err', 'flux', 'flux_err':
                spectra[order][key][i] = s[order][key]

        hdulist.close()

    retval = { }
    for order in spectra:
        weights = spectra[order]['livetime'] / spectra[order]['livetime'].sum()
        weights = weights[:,np.newaxis]

        retval[order] = { }
        retval[order]['lambda'] = wav

        retval[order]['flux'] = (weights * spectra[order]['flux']).sum(axis=0)
        retval[order]['flux_err'] = np.sqrt((weights**2 * spectra[order]['flux_err']**2).sum(axis=0))

        retval[order]['rate'] = (weights * spectra[order]['rate']).sum(axis=0)
        retval[order]['rate_err'] = np.sqrt((weights**2 * spectra[order]['rate_err']**2).sum(axis=0))

    return retval

def bins_to_energy(bin_lo, bin_hi):
    return 12.39854 / (0.5 * (bin_lo + bin_hi))
