#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Create precip type maps from WRF output

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

# Define function to determine precip type
def ptype(frac,reflectivity,temp):
    # List containing fraction of frozen precip
    type_pct = frac
        
    # Set reflectivity level for rearrangement
    refl = reflectivity
        
    # Empty lists to contain reflectivity data
    rain_precip = []
    zr_precip = []
    slt_precip = []
    snow_precip = []
    
    # Loop through all lons
    for j in range(len(refl)):
        # Empty lists to contain individual column data
        cur_col_rain = []
        cur_col_zr = []
        cur_col_slt = []
        cur_col_snow = []
        # Loop through all lats
        for k in range(len(refl[0])):
            # Check for certain fraction ranges
            if type_pct[j,k] >= 0.9:
                cur_col_rain.append(0.)
                cur_col_zr.append(0.)
                cur_col_slt.append(0.)
                cur_col_snow.append(refl[j,k])
            elif (0.1 < type_pct[j,k] < 0.9):
                cur_col_rain.append(0.)
                cur_col_zr.append(0.)
                cur_col_slt.append(refl[j,k])
                cur_col_snow.append(0.)
            elif type_pct[j,k] < 0.1 and temp[j,k] < 273.15:
                cur_col_rain.append(0.)
                cur_col_zr.append(refl[j,k])
                cur_col_slt.append(0.)
                cur_col_snow.append(0.)
            elif type_pct[j,k] < 0.1 and temp[j,k] >= 273.15:
                cur_col_rain.append(refl[j,k])
                cur_col_zr.append(0.)
                cur_col_slt.append(0.)
                cur_col_snow.append(0.)
        # Append each column value to new lists for plotting
        snow_precip.append(cur_col_snow)
        slt_precip.append(cur_col_slt)
        zr_precip.append(cur_col_zr)
        rain_precip.append(cur_col_rain) 
        
    return snow_precip, slt_precip, zr_precip, rain_precip

# Open the NetCDF file
inDir = "../WRF/run/"
inFile = inDir + "wrfout_d02_2022-02-14_00:00:00"
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
os.system("rm ./images/ptype/*")

#### Loop through each timestep and create a new plot for each one
for i in range(len(dbz)):
    
    # The first timestep is only good for setting the initialization time for 
    # the plot title. Afterwards, only plot the subsequent timesteps.
    if i == 0:
        init = str(times[i])[11:13] + 'z ' + str(times[i])[0:10]
        
    else:
        # Change the DPI of the resulting figure. Higher DPI drastically improves the
        # look of the text rendering.
        plt.rcParams['savefig.dpi'] = 255
        
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
                     transform=crs.PlateCarree(), zorder=1, alpha=0.5)
        
        # Use the line contours to place contour labels.
        ax.clabel(
            line_c,  # Typically best results when labelling line contours.
            colors=['black'],
            manual=False,  # Automatic placement vs manual placement.
            inline=True,  # Cut the line where the label will be placed.
            fmt=' {:.0f} '.format,  # Labels as integers, with some extra space.
            fontsize=5,
        )
        
        ## Loop through the grid to determine the ptype at each grid location
        snowPrecip, sltPrecip, zrPrecip, rainPrecip = ptype(sr[i],dbz[i][0],t2[i])
        
        # Make filled contour for the dBZ values for rain, ZR, IP, and snow
        # The colormaps are custom, so feel free to change them as you choose
        dbzRng = np.arange(0,76,5)
        rainMesh = plt.contourf(to_np(lons), to_np(lats), to_np(rainPrecip), dbzRng,
                     transform=crs.PlateCarree(), zorder=0, cmap=cm.dbz_colormap, extend='max')
        
        zrMesh = plt.contourf(to_np(lons), to_np(lats), to_np(zrPrecip), dbzRng,
                     transform=crs.PlateCarree(), zorder=0, cmap=cm.zr_colormap, extend='max')
        
        sltMesh = plt.contourf(to_np(lons), to_np(lats), to_np(sltPrecip), dbzRng,
                     transform=crs.PlateCarree(), zorder=0, cmap=cm.sleet_colormap, extend='max')
        
        snowMesh = plt.contourf(to_np(lons), to_np(lats), to_np(snowPrecip), dbzRng,
                     transform=crs.PlateCarree(), zorder=0, cmap=cm.snow_colormap, extend='max')
    
        # Add a color bar
        cax = plt.colorbar(rainMesh, shrink=.65, ax=ax, pad=-0.06)
        cax2 = plt.colorbar(zrMesh, shrink=.65, ax=ax, pad=-0.05)
        cax3 = plt.colorbar(sltMesh, shrink=.65, ax=ax, pad=-0.03)
        cax4 = plt.colorbar(snowMesh, shrink=.65, ax=ax)
        cax.ax.tick_params(labelsize=6)
        cax2.ax.tick_params(labelsize=6)
        cax3.ax.tick_params(labelsize=6)
        cax4.ax.tick_params(labelsize=6)
        cax.ax.set_title('dBZ', fontsize=6)
        cax2.ax.set_title('ZR', fontsize=6)
        cax3.ax.set_title('SLEET', fontsize=6)
        cax4.ax.set_title('SNOW', fontsize=6)
    
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
        title_line1 = str('WRF-3km Simulated P-Type Reflectivity (fill) and MSLP (contour, hPa)')
        title_line2 = str('Initialized: '+init+'      Time Step: '+ step +'      Valid at: '+ hour + ' ' + date)
        plt.title(title_line1+'\n'+title_line2, fontsize=8)
        
        # Add additional info to the plot (creator info, max/min values, etc.)
        ax.text(0,-0.01,'Matthew Toadvine\nTwitter: @toadwx', horizontalalignment='left',
                verticalalignment='top', transform=ax.transAxes, fontsize=5)
        #ax.text(1,-0.01,'Maximum 2m Temperature: {:5.1f}\N{DEGREE SIGN}F\nMinimum Pressure: {:6.1f} hPa'.format(maxTemp,minMSLP), 
        #    horizontalalignment='right', verticalalignment='top', transform=ax.transAxes)
        
        # Save the figures
        plt.savefig('./images/ptype/step_'+step+'.png', format='png',bbox_inches='tight')
        #plt.show()
        plt.close()
        
        # Append the image to the list frames for GIF processing
        image = './images/ptype/step_'+step+'.png'
        frames.append(Image.open(image))

# Perform automatic GIF processing and save the file
frame_one = frames[0]
frame_one.save("./images/ptype/ptype.gif", format="GIF", append_images=frames,
           save_all=True, duration=500, loop=0)
