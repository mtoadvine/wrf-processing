#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Create 2 meter Temperature plots from WRF output.

@author: Matthew Toadvine
@date:   01/14/2022
"""

# Import required libraries
import sys, os
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import cartopy.crs as crs
import cartopy.feature as cfeature
import metpy.calc as mpcalc
import metpy.units as units
from xarray.backends import NetCDF4DataStore
import xarray as xr
from math import *
from wrf import (to_np, getvar, smooth2d, get_cartopy, cartopy_xlim,
                 cartopy_ylim, latlon_coords, ALL_TIMES, ll_to_xy)

# Open the NetCDF file
inDir = "../WRF/run/"
inFile = inDir + "wrfout_d01_2022-02-07_00:00:00"
ncfile = Dataset(inFile)

# Set desired sounding Lat/Lon 
lat = 35.1919
lon = -80.6827

# Check to see if the Lat/Lon point is within the domain
# If invalid, quit with error message
try:
    x, y = ll_to_xy(ncfile,lat,lon)
except:
    sys.exit('Lat/Lon point is not within the domain.')

# Get the variable of interest and timesteps of interest
temp2m = getvar(ncfile, "T2", timeidx=ALL_TIMES, method='cat')
slp = getvar(ncfile, "slp", timeidx=ALL_TIMES, method='cat')
times = list(temp2m.coords['Time'].values)

# Convert the temperature array from Kelvin to deg C/F 
# This requires converting the xarray dataset to quantify units first
# Conversion optional - be sure to edit plot texts accordingly
temp2m = temp2m.metpy.quantify()
#temp2m = temp2m.metpy.convert_units('degC')
temp2m = temp2m.metpy.convert_units('degF')

# These lines identify the max/min temps within the dataset to provide an
# appropriate range for the colortable when plotting
rngMin = floor((np.amin((temp2m.squeeze()).metpy.dequantify()))/10)*10
rngMax = ceil((np.amax((temp2m.squeeze()).metpy.dequantify()))/10)*10

# Smooth the sea level pressure since it tends to be noisy near the
# mountains
smooth_slp = smooth2d(slp, 3, cenweight=4)

# Get the latitude and longitude points
lats, lons = latlon_coords(slp)

# Get the cartopy mapping object
cart_proj = get_cartopy(slp)

# Create a list to hold all images files for GIF processing
# and remove all old images from the directory
frames = []
os.system("rm ./images/temp2m/*")

#### Loop through each timestep and create a new plot for each timestep
for i in range(len(temp2m)):
    
    # The first timestep is only good for setting the initialization time for 
    # the plot title. Afterwards, only plot the subsequent timesteps.
    if i == 0:
        init = str(times[i])[11:13] + 'z ' + str(times[i])[0:10]
        
    else:
        # Change the DPI of the resulting figure. Higher DPI drastically improves the
        # look of the text rendering.
        plt.rcParams['savefig.dpi'] = 255        
        
        # Create a new figure
        fig = plt.figure(figsize=(12,10))
        
        # Set the GeoAxes to the projection used by WRF
        ax = plt.axes(projection=cart_proj)
    
        # Add the states and coastlines
        ax.add_feature(cfeature.STATES, edgecolor='black', linewidth=0.5)
        ax.add_feature(cfeature.COASTLINE, linewidth=0.8)
        ax.coastlines('50m', linewidth=0.8)
        
        # Make filled contour for the temperature values using a defined range
        # based on the max/min values
        tempRng = np.arange(rngMin,rngMax,2)
        tempMesh = plt.contourf(to_np(lons), to_np(lats), to_np(temp2m[i]), tempRng,
                     transform=crs.PlateCarree(), zorder=0, cmap=get_cmap("jet"), extend='both')
        frzLine = plt.contour(to_np(lons), to_np(lats), to_np(temp2m[i]), np.arange(32,33),
                               transform=crs.PlateCarree(), zorder=0, colors="black", linestyles='dashed')
    
        # Add a color bar
        cbar = fig.colorbar(tempMesh, shrink=.98)
        cbar.set_label('\N{DEGREE SIGN}F')
    
        # Set the map boundaries based on data extent
        ax.set_xlim(cartopy_xlim(smooth_slp))
        ax.set_ylim(cartopy_ylim(smooth_slp))
    
        # Add the gridlines
        ax.gridlines(color="black", linestyle="dotted")
        
        # Find minimum MSLP value for display
        #minTemp = np.amin(temp2m[i].squeeze())
        #minTemp = minTemp.values
    
        # Find location-specific (lat/lon) temperature for display
        selectTemp = temp2m[i][y][x].values
        
        
        #### Add plot titles and basic model run information
        
        # Define timestep hour and date for title    
        date = str(times[i])[0:10]
        hour = str(times[i])[11:16] + 'z'
        if i < 10:
            step = '0' + str(i)
        else:
            step = str(i)
        
        # Create two title lines
        title_line1 = str('WRF-9km 2m Temperature')
        title_line2 = str('Initialized: '+init+'      Time Step: '+ step +'      Valid at: '+ hour + ' ' + date)
        plt.title(title_line1+'\n'+title_line2)
        
        # Add additional info to the plot (creator info, max/min values, etc.)
        # These go underneath the main plot
        ax.text(0,-0.01,'Matthew Toadvine\nTwitter: @toadwx', horizontalalignment='left',
                verticalalignment='top', transform=ax.transAxes)
        ax.text(1,-0.01,'Charlotte 2m Temperature: {:5.1f}\N{DEGREE SIGN}F'.format(selectTemp), 
            horizontalalignment='right', verticalalignment='top', transform=ax.transAxes)
        
        # Save the figures
        plt.savefig('./images/temp2m/step_'+step+'.png', format='png',bbox_inches='tight')
        #plt.show()
        plt.close()
