import os
import sys

import pylab
import numpy
from mpl_toolkits.basemap import Basemap


map = Basemap(projection='ortho',lat_0=0,lon_0=-70,
              resolution='l',area_thresh=1000.)
#map.drawcoastlines()
#map.drawcountries()
map.bluemarble()
#map.fillcontinents(color='coral')
#map.drawmapboundary()
#map.drawmeridians(numpy.arange(0,360,30))
#map.drawparallels(numpy.arange(-90,90,30))

pylab.savefig("globe.jpg", fmt="jpg", dpi=200)
