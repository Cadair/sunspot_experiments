import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from skimage.filters import median
from skimage.measure import label, regionprops
from skimage.morphology import opening, square, watershed, white_tophat, disk
from skimage.util import invert

import sunpy.map

map1 = sunpy.map.Map("./data/hmi_ic_45s_2014_02_13_00_01_30_tai_continuum.fits")
map1 = map1.resample((1024, 1024)*u.pix)

data = invert(map1.data)

circle_radius = 40 * u.arcsec
median_box = 10 * u.arcsec

c_pix = int((circle_radius / map1.scale[0]).to_value(u.pix))
circle = disk(c_pix/2)

m_pix = int((median_box / map1.scale[0]).to_value(u.pix))
med = median(data, square(m_pix), behavior="ndimage")
th = white_tophat(med, circle)

finite = th
finite[np.isnan(th)] = 0

segmentation = finite > np.percentile(finite, 99.8)
labelled = label(segmentation)
regions = regionprops(labelled)

plt.figure()
ax = plt.subplot(projection=map1)
map1.plot()
ax.contour(segmentation, levels=0)
