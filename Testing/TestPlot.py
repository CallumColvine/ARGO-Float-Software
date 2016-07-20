# Plotting Imports
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
# Calculations imports
import numpy as np

# setup polyconic basemap
# by specifying lat/lon corners and central point.
# area_thresh=1000 means don't plot coastline features less
# than 1000 km^2 in area.
m = Basemap(llcrnrlon=-160,llcrnrlat=43,urcrnrlon=-100,urcrnrlat=57,\
            resolution='l',area_thresh=1000.,projection='poly',\
            lat_0=50,lon_0=-140)
m.drawcoastlines()
m.fillcontinents(color='darksage',lake_color='royalblue')
# draw parallels and meridians.
m.drawparallels(np.arange(-80.,81.,20.))
m.drawmeridians(np.arange(-180.,181.,20.))
m.drawmapboundary(fill_color='royalblue')
plt.title("Circulation Data on the Southern Alaskan Coast")
for x in xrange(50, 65):
    plt.plot(-155, x)
    
plt.show()