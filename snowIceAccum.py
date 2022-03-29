#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Create snowfall maps from WRF output

@author: Matthew Toadvine
@date:   01/21/2022
"""

# Import required libraries
from netCDF4 import Dataset
import numpy as np
from math import *
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import cartopy.crs as crs
import cartopy.feature as cfeature
import colormap as cm
from wrf import (to_np, getvar, smooth2d, get_cartopy, cartopy_xlim,
                 cartopy_ylim, latlon_coords, ALL_TIMES)


# Open the NetCDF file
inDir = "../WRF/run/"
inFile = inDir + "wrfout_d01_2022-01-29_00:00:00"
ncFile = Dataset(inFile)


# Get the variable of interest and timesteps of interest
snownc = getvar(ncFile, "SNOWNC", timeidx=ALL_TIMES, method='cat')
slp = getvar(ncFile, "slp", timeidx=ALL_TIMES, method='cat')
times = list(snownc.coords['Time'].values)

# Create an empty NP array of zeroes the size of 

# These lines identify the max/min temps within the dataset to provide an
# appropriate range for the colortable when plotting
rngMin = floor(np.amin(snownc.values))
rngMax = ceil(np.amax(snownc.values))

# Smooth the sea level pressure since it tends to be noisy near the
# mountains
smooth_slp = smooth2d(slp, 3, cenweight=4)

# Get the latitude and longitude points
lats, lons = latlon_coords(slp)

# Get the cartopy mapping object
cart_proj = get_cartopy(slp)
       
        

#### Loop through each timestep and create a new plot for each one
for i in range(len(snownc)):
    
    # The first timestep is only good for setting the initialization time for 
    # the plot title. Afterwards, only plot the subsequent timesteps.
    if i == 0:
        init = str(times[i])[11:13] + 'z ' + str(times[i])[0:10]
        
    else:
        # Change the DPI of the resulting figure. Higher DPI drastically improves the
        # look of the text rendering.
        plt.rcParams['savefig.dpi'] = 255        
        
        # Create a figure
        fig = plt.figure(figsize=(20,12))
        
        # Set the GeoAxes to the projection used by WRF
        ax = plt.axes(projection=cart_proj)
    
        # Add the states and coastlines
        ax.add_feature(cfeature.STATES, edgecolor='black', linewidth=0.5)
        ax.add_feature(cfeature.COASTLINE, linewidth=0.8)
        ax.coastlines('50m', linewidth=0.8)
        
        # Make filled contour for the dBZ values for rain, mix, and snow
        snowLevels = np.arange(0,5,0.25)
        snowMesh = plt.contourf(to_np(lons), to_np(lats), to_np(snownc[i]), snowLevels,
                     transform=crs.PlateCarree(), zorder=0, cmap=get_cmap("BuPu"), extend='both')
    
        # Add a color bar
        cax = plt.colorbar(snowMesh, shrink=.65)
        cax.set_label('inches')
    
        # Set the map bounds
        ax.set_xlim(cartopy_xlim(smooth_slp))
        ax.set_ylim(cartopy_ylim(smooth_slp))
    
        # Add the gridlines
        #ax.gridlines(color="black", linestyle="dotted")
        
        # Find minimum MSLP value for display
        #minMSLP = np.amin(slp[i].squeeze())
        #minMSLP = minMSLP.values
    
        # Find temperature for display
        #maxTemp = np.amax(temp2m[i].squeeze())
        #maxTemp = maxTemp.values
        
        
        #### Add plot titles and basic model run information
        
        # Define timestep hour and date for title    
        date = str(times[i])[0:10]
        hour = str(times[i])[11:16] + 'z'
        if i < 10:
            step = '0' + str(i)
        else:
            step = str(i)
        
        # Create two title lines
        title_line1 = str('WRF-9km Simulated Accumulated Snowfall (fill)')
        title_line2 = str('Initialized: '+init+'      Time Step: '+ step +'      Valid at: '+ hour + ' ' + date)
        plt.title(title_line1+'\n'+title_line2)
        
        # Add additional info to the plot (creator info, max/min values, etc.)
        ax.text(0,-0.01,'Matthew Toadvine\nTwitter: @toadwx', horizontalalignment='left',
                verticalalignment='top', transform=ax.transAxes)
        #ax.text(1,-0.01,'Maximum 2m Temperature: {:5.1f}\N{DEGREE SIGN}F\nMinimum Pressure: {:6.1f} hPa'.format(maxTemp,minMSLP), 
        #    horizontalalignment='right', verticalalignment='top', transform=ax.transAxes)
        
        # Save the figures
        plt.savefig('./images/snow/step_'+step+'.png', format='png',bbox_inches='tight')
        #plt.show()
        plt.close()
