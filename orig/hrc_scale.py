import numpy as np
import os, errno

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
