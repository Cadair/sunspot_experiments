from functools import partial

import astropy.units as u
import matplotlib.pyplot as plt

import sunpy.map

from stara import stara

maps = sunpy.map.Map("./data/*720*")
maps = [m.resample((1024, 1024) * u.pix) for m in maps]

segs = list(map(partial(stara, limb_filter=10*u.percent), maps))

for smap, seg in zip(maps, segs):
    plt.figure()
    ax = plt.subplot(projection=smap)
    smap.plot()
    ax.contour(seg, levels=0)
    plt.colorbar()

plt.show()
