import numpy as np

def combine_fluxes(data, inst):
    iip = np.where((data['inst']==inst) & (data['order']=='pos'))[0]
    iin = np.where((data['inst']==inst) & (data['order']=='neg'))[0]

    if np.where(data['days'][iip] != data['days'][iin])[0].size:
        raise ValueError("uh oh")

    flux = data['flux'][iip] + data['flux'][iin]
    flux_neg = data['flux'][iin]
    flux_pos = data['flux'][iip]

    # "fits" RDB flies have flux_min, flux_max columns
    if 'flux_max' in data:
        flux_err = np.sqrt(
            (0.5*(data['flux_max'][iip]-data['flux_min'][iip]))**2 +
            (0.5*(data['flux_max'][iin]-data['flux_min'][iin]))**2
        )
        flux_neg_err = 0.5*(data['flux_max'][iin]-data['flux_min'][iin])
        flux_pos_err = 0.5*(data['flux_max'][iip]-data['flux_min'][iip])

    # "obsflux" files have err column
    else:
        flux_err = np.sqrt( data['err'][iip]**2 + data['err'][iin]**2 )
        flux_neg_err = data['err'][iin]
        flux_pos_err = data['err'][iip]

    return { 'obsid' : data['obsid'][iip],
             'days'  : data['days'][iip],
             'flux'  : flux,
             'flux_err'   : flux_err,
             'flux_neg' : flux_neg,
             'flux_neg_err' : flux_neg_err,
             'flux_pos' : flux_pos,
             'flux_pos_err' : flux_pos_err,
    }
