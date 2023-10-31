from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np

# coordinates of Hannover, Germany
lat = 52.3758916
lon = 9.732010400000038

for name in ["topo_no_lines", "topo", "plain"]:
    plt.clf()

    # get world map
    map = Basemap(projection='ortho',lat_0=lat,lon_0=lon,resolution='l')
    map.drawcoastlines(linewidth=0.25)
    map.drawcountries(linewidth=0.25)
    if name != "topo_no_lines":
        map.drawparallels(range(-90,91,30))
        map.drawmeridians(range(-180,181,30))

    # draw topoligical features or plain map (water and land)
    if "topo" in name:
        map.etopo()
        if name == "topo_no_lines":
            # draw x, y and z axis

            # draw a point at the north pole
            x, y = map(0, 90.)
            plt.plot(x, y, 'o', color='red', markersize=3)

            # draw a point for end of z axis
            x1 = x
            y1 = y + 1100000
            # plot arrow from x,y to x1, y1
            plt.arrow(x, y, x1-x, y1-y, head_width=100000, head_length=100000, color="red")
            plt.text(x1+150000, y1-100000, "z", color="red", fontsize=10)

            # draw a point at lon = 0, lat = 0
            x, y = map(0, 0)
            plt.plot(x, y, 'o', color='red', markersize=3)

            # draw a point for end of x axis
            x2 = x -  200000
            y2 = y - 1300000
            # plot arrow from x,y to x2, y2
            plt.arrow(x, y, x2-x, y2-y, head_width=100000, head_length=100000, color="red")
            plt.text(x2-350000, y2+200000, "x", color="red", fontsize=10)

            # draw a point at lon = 90, lat = 0
            x, y = map(90., 0)
            plt.plot(x, y, 'o', color='red', markersize=3)

            # draw a point for end of y axis
            x3 = x + 1300000
            y3 = y - 250000
            # plot arrow from x,y to x3, y3
            plt.arrow(x, y, x3-x, y3-y, head_width=100000, head_length=100000, color="red")
            plt.text(x3-450000,y3+250000, "y", color="red", fontsize=8)

            # get x-lim and y-lim
            xlim = plt.gca().get_xlim()
            ylim = plt.gca().get_ylim()

            plt.gca().set_xlim((xlim[0]*1.1, xlim[1]*1.15))
            plt.gca().set_ylim((ylim[0]*1.1, ylim[1]*1.1))

    elif name == "plain":
        map.fillcontinents(color='tan',lake_color='lightgray')
        map.drawmapboundary(fill_color='lightgray')
    
    plt.savefig("globe_" + name + ".pdf")