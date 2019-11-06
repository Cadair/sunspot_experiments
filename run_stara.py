from functools import partial

import astropy.units as u
import matplotlib.pyplot as plt
from astropy.time import Time
from astropy.table import QTable, vstack
from skimage.measure import regionprops, label, regionprops_table
from skimage.color import label2rgb

import sunpy.map

from stara import stara

maps = sunpy.map.Map("./data/*720*")
maps = [m.resample((1024, 1024) * u.pix) for m in maps]

segs = list(map(partial(stara, limb_filter=10 * u.percent), maps))

def get_regions(segmentation, smap):
    labelled = label(segmentation)
    if labelled.max() == 0:
        return QTable()

    regions = regionprops_table(labelled, smap.data,
                                properties=["label",
                                            "centroid",
                                            "area",
                                            "min_intensity"])

    regions['obstime'] = Time([smap.date] * regions['label'].size)
    regions['center_coord'] = smap.pixel_to_world(regions['centroid-0'] * u.pix,
                                                  regions['centroid-1'] * u.pix).heliographic_stonyhurst

    return QTable(regions)


print(list(map(get_regions, segs, maps)))


# for smap, seg in zip(maps, segs):
#     plt.figure()
#     ax = plt.subplot(projection=smap)
#     smap.plot()
#     ax.contour(seg, levels=0)
#     plt.colorbar()

plt.show()
