import astropy.units as u

from sunpy.net import Fido
from sunpy.net import attrs as a
from sunpy.time import parse_time
from sunpy.net.attr import or_


times = map(parse_time, ["2014-02-13", "2019-08-01"])
times = list(map(lambda t: a.Time(t, t+1*u.day, t), times))
times = or_(*times)

results = Fido.search(a.Instrument("HMI"),
                      a.vso.Physobs("intensity"),
                      times)

print(results)

files = Fido.fetch(results, path="./data/{file}")
files = sorted(files)
print(files)
