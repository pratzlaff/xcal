import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc, rcParams
rc('text', usetex=True)
rcParams.update({'font.size': 14})

wav, data, err, model = np.loadtxt('/data/legs/rpete/flight/xcal/mkn421_hrcs_letg_coadded/foo.txt', unpack=True)

f, axarr = plt.subplots(2, sharex=True, figsize=(8.5, 11))

x = wav
y = data
yerr = err

axarr[0].plot(x, y, 'k.')
axarr[0].errorbar(x, y, yerr, ecolor='k', fmt=None)
axarr[0].plot(x, model, 'r-')
axarr[0].set_ylabel(r'Negative 1st Order Data and Model')
axarr[0].set_title('Mkn 421 HRC-S/LETG, coadded')
axarr[0].set_yscale('log')

x = wav
y = data/model
yerr = err/model

axarr[1].plot(x, y, 'ko')
axarr[1].errorbar(x, y, yerr, ecolor='k', fmt=None)
axarr[1].set_xlim(42, 44)
axarr[1].set_ylim(0.4, 1.3)
axarr[1].set_xlabel(r'$\lambda$ (\AA)')
axarr[1].set_ylabel(r'Negative 1st Order Data/Model Ratio')

f.tight_layout()
f.subplots_adjust(hspace=0)

plt.show()
