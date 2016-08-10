# # Plotting Imports
# from mpl_toolkits.basemap import Basemap
# import matplotlib.pyplot as plt
# # Calculations imports
# import numpy as np

# # setup polyconic basemap
# # by specifying lat/lon corners and central point.
# # area_thresh=1000 means don't plot coastline features less
# # than 1000 km^2 in area.

# m = Basemap(llcrnrlon=-160,llcrnrlat=43,urcrnrlon=-100,urcrnrlat=57,\
#             resolution='l',area_thresh=1000.,projection='poly',\
#             lat_0=50,lon_0=-140)

# # m = Basemap(llcrnrlon=-160,llcrnrlat=43,urcrnrlon=-100,urcrnrlat=57,\
# #             resolution='l',area_thresh=1000.,projection='poly',\
# #             lat_0=50,lon_0=-140)


# m.drawcoastlines()
# m.fillcontinents(color='darksage',lake_color='royalblue')
# # draw parallels and meridians.
# m.drawparallels(np.arange(-80.,81.,20.))
# m.drawmeridians(np.arange(-180.,181.,20.))
# m.drawmapboundary(fill_color='royalblue')
# plt.title("Circulation Data on the Southern Alaskan Coast")

# # for x in xrange(50, 65):
# #     plt.plot(-155, x)
    
# plt.show()




from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
# set up orthographic map projection with
# perspective of satellite looking down at 50N, 100W.
# use low resolution coastlines.
map = Basemap(projection='ortho',lat_0=45,lon_0=-100,resolution='l')
# draw coastlines, country boundaries, fill continents.
map.drawcoastlines(linewidth=0.25)
map.drawcountries(linewidth=0.25)
map.fillcontinents(color='coral',lake_color='aqua')
# draw the edge of the map projection region (the projection limb)
map.drawmapboundary(fill_color='aqua')
# draw lat/lon grid lines every 30 degrees.
map.drawmeridians(np.arange(0,360,30))
map.drawparallels(np.arange(-90,90,30))
# make up some data on a regular lat/lon grid.
nlats = 73; nlons = 145; delta = 2.*np.pi/(nlons-1)
lats = (0.5*np.pi-delta*np.indices((nlats,nlons))[0,:,:])
lons = (delta*np.indices((nlats,nlons))[1,:,:])
wave = 0.75*(np.sin(2.*lats)**8*np.cos(4.*lons))
mean = 0.5*np.cos(2.*lats)*((np.sin(2.*lats))**2 + 2.)

print "len lats ", len(lats), " len lons ", len(lons)
print "len lats[0] ", len(lats[0]), " len lons[0] ", len(lons[0])
print "len wave/mean ", len(wave), "len wave[0]", len(wave[0])

print "lats", lats
print "lons", lons

# print "Wave: ", wave
# print "Mean: ", mean
# print "Wave + Mean: ", wave + mean
# compute native map projection coordinates of lat/lon grid.
x, y = map(lons*180./np.pi, lats*180./np.pi)
# contour data over the map.
cs = map.contour(x,y,wave+mean,15,linewidths=1.5)
plt.title('contour lines over filled continent background')
plt.show()