#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Create radar reflectivity maps from WRF output

@author: Matthew Toadvine
@date:   01/14/2022
"""

# Import required libraries
from netCDF4 import Dataset
import os
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as crs
import cartopy.feature as cfeature
import colormap as cm
from PIL import Image
from wrf import (to_np, getvar, smooth2d, get_cartopy, cartopy_xlim,
                 cartopy_ylim, latlon_coords, ALL_TIMES)

# Open the NetCDF file
inDir = "../WRF/run/"
inFile = inDir + "wrfout_d01_2022-03-23_12:00:00"
ncFile = Dataset(inFile)

# Get the variable of interest and timesteps of interest
#t2 = getvar(ncfile, 'T2', timeidx=ALL_TIMES, method='cat')
dbz = getvar(ncFile, "dbz", timeidx=ALL_TIMES, method='cat')
slp = getvar(ncFile, "slp", timeidx=ALL_TIMES, method='cat')
sr = getvar(ncFile, "SR", timeidx=ALL_TIMES, method='cat')
t2 = getvar(ncFile, "T2", timeidx=ALL_TIMES, method='cat')
times = list(dbz.coords['Time'].values)


# These lines identify the max/min temps within the dataset to provide an
# appropriate range for the colortable when plotting
#rngMin = floor((np.amin((temp2m.squeeze()).metpy.dequantify()))/10)*10
#rngMax = ceil((np.amax((temp2m.squeeze()).metpy.dequantify()))/10)*10

# Smooth the sea level pressure since it tends to be noisy near the
# mountains
smooth_slp = smooth2d(slp, 3, cenweight=4)

# Get the latitude and longitude points
lats, lons = latlon_coords(slp)

# Get the cartopy mapping object
cart_proj = get_cartopy(slp)
       
# Define a list to contain all of the image files for automatic 
# GIF processing and remove all old images from the directory
frames = []       
os.system("rm ./images/refl/*")

#### Loop through each timestep and create a new plot for each one
for i in range(len(dbz)):
    
    # The first timestep is only good for setting the initialization time for 
    # the plot title. Afterwards, only plot the subsequent timesteps.
    if i == 0:
        init = str(times[i])[11:13] + 'z ' + str(times[i])[0:10]
        
    else:
        # Change the DPI of the resulting figure. Higher DPI drastically improves the
        # look of the text rendering.
        plt.rcParams['savefig.dpi'] = 300
        
        # Create a figure
        fig = plt.figure(figsize=(10,6))
        
        # Set the GeoAxes to the projection used by WRF
        ax = plt.axes(projection=cart_proj)
    
        # Add the states and coastlines
        ax.add_feature(cfeature.STATES, edgecolor='black', linewidth=0.5)
        ax.add_feature(cfeature.COASTLINE, linewidth=0.8)
        ax.coastlines('50m', linewidth=0.8)
    
        # Make the contour outlines for smoothed SLP values
        slpRng = np.arange(950,1051,4)
        line_c = plt.contour(to_np(lons), to_np(lats), to_np(smooth_slp[i]), slpRng, colors="black",
                     transform=crs.PlateCarree(), zorder=1, alpha=0.3, linewidths=0.3)
        
        # Use the line contours to place contour labels.
        ax.clabel(
            line_c,  # Typically best results when labelling line contours.
            colors=['black'],
            manual=False,  # Automatic placement vs manual placement.
            inline=True,  # Cut the line where the label will be placed.
            fmt=' {:.0f} '.format,  # Labels as integers, with some extra space.
            fontsize=5,
        )
        
        # Make filled contour for the dBZ values for rain, ZR, IP, and snow
        # The colormaps are custom, so feel free to change them as you choose
        dbzRng = np.arange(0,76,5)
        rainMesh = plt.contourf(to_np(lons), to_np(lats), to_np(dbz[i][0]), dbzRng,
                     transform=crs.PlateCarree(), zorder=0, cmap=cm.dbz_colormap, extend='max')
    
        # Add a color bar
        cax = plt.colorbar(rainMesh, shrink=.65, ax=ax)
        cax.ax.tick_params(labelsize=6)
        cax.ax.set_title('dBZ', fontsize=6, y=1.05)
        #cax.set_label('dBZ', loc='bottom')
    
        # Set the map bounds
        ax.set_xlim(cartopy_xlim(smooth_slp))
        ax.set_ylim(cartopy_ylim(smooth_slp))
    
        # Add the gridlines
        #ax.gridlines(color="black", linestyle="dotted")
        
        
        #### Add plot titles and basic model run information
        
        # Define timestep hour and date for title    
        date = str(times[i])[0:10]
        hour = str(times[i])[11:16] + 'z'
        if i < 10:
            step = '0' + str(i)
        else:
            step = str(i)
        
        # Create two title lines
        title_line1 = str('WRF-3km Simulated Radar Reflectivity (fill) and MSLP (contour, hPa)')
        title_line2 = str('Initialized: '+init+'      Time Step: '+ step +'      Valid at: '+ hour + ' ' + date)
        plt.title(title_line1+'\n'+title_line2, fontsize=8)
        
        # Add additional info to the plot (creator info, max/min values, etc.)
        ax.text(0,-0.01,'Matthew Toadvine\nTwitter: @toadwx', horizontalalignment='left',
                verticalalignment='top', transform=ax.transAxes, fontsize=5)
        #ax.text(1,-0.01,'Maximum 2m Temperature: {:5.1f}\N{DEGREE SIGN}F\nMinimum Pressure: {:6.1f} hPa'.format(maxTemp,minMSLP), 
        #    horizontalalignment='right', verticalalignment='top', transform=ax.transAxes)
        
        # Save the figures
        plt.savefig('./images/refl/step_'+step+'.png', format='png',bbox_inches='tight')
        #plt.show()
        plt.close()
        
        # Append the image to the list frames for GIF processing
        image = './images/refl/step_'+step+'.png'
        frames.append(Image.open(image))

# Perform automatic GIF processing and save the file
frame_one = frames[0]
frame_one.save("./images/refl/refl.gif", format="GIF", append_images=frames,
           save_all=True, duration=500, loop=0)
