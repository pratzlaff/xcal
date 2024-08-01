#! /bin/bash

crop="15 2 140 4 wcs"

ds9 -geometry 900x1200 \
  -colorbar no \
  -zoom 2 \
  -scale mode 99.5 \
  -grid yes \
  -grid type publication \
  -grid labels yes \
  -grid title gap 10 \
  -grid title def no \
  -grid labels def1 no \
  -grid labels def2 no \
  -cmap Heat \
  data/mkn421/16429/osort_img_osip_leg_all.fits \
  data/mkn421/16429/osort_img_osip_leg_123.fits \
  data/mkn421/16429/osort_img_noosip_leg_all.fits \
  data/mkn421/16429/osort_img_noosip_leg_123.fits \
  -tile mode grid \
  -frame 1 -crop $crop -grid title text 'OSIP: All Orders' -grid labels text2 Order \
  -frame 2 -crop $crop -grid title text 'OSIP: Orders 1-3' \
  -frame 3 -crop $crop -grid title text 'No OSIP: All Orders' -grid labels text1 TG_MLAM -grid labels text2 Order \
  -frame 4 -crop $crop -grid title text 'No OSIP: Orders 1-3' -grid labels text1 TG_MLAM \
  -print destination file \
  -print filename ds9.ps \
  -print

