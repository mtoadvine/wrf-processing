#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Algorithm and plotting of Austin Slip Index (ASI) from WRF output

@author: Matthew Toadvine
@date:   01/15/2022
"""

# Import required libraries
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import cartopy.crs as crs
import cartopy.feature as cfeature
import colormap as cm
from wrf import (to_np, getvar, get_cartopy, cartopy_xlim,
                 cartopy_ylim, latlon_coords, ALL_TIMES)

# Define function to calculate the ASI
def slipIndex(time, temp, frac, reflectivity, uwind, vwind):
    # Convert datetime object to integer of hour
    hr = int(str(time)[11:13])
        
    # Rename variables for simplicity
    refl = reflectivity
    speed = np.sqrt(uwind**2 + vwind**2)
    
    # Define max values for normalization prior to weights
    refMax = 75.
    spdMax = 50.
    
    # Define static weights
    fracWt = 0.6
    windWt = 0.2
    reflWt = 0.2
        
    # Pre-allocate the array to hold our ASI values
    slip = np.zeros(shape=(len(refl),len(refl[0])))
    
    # Loop through all lons
    for j in range(len(refl)):
        # Loop through all lats
        for k in range(len(refl[0])):
            # Check for certain conditions and adjust certain weights as necessary
            if temp[j,k] < 273.15:
                if hr > 12:
                    timeWt = 0.05
                    fracWt = 0.55
                    weight = timeWt + fracWt + ((refl[j,k]/refMax)*reflWt) + ((speed[j,k]/spdMax)*windWt)
                    slip[j,k] = weight
                else:
                    timeWt = 0.
                    fracWt = 0.6
                    weight = timeWt + fracWt + ((refl[j,k]/refMax)*reflWt) + ((speed[j,k]/spdMax)*windWt)
                    slip[j,k] = weight
            else:
                if hr > 12:
                    timeWt = 0.05
                    fracWt = 0.55
                    weight = 0.1 * (timeWt + fracWt + ((refl[j,k]/refMax)*reflWt) + ((speed[j,k]/spdMax)*windWt))
                    slip[j,k] = weight
                else:
                    timeWt = 0.
                    fracWt = 0.6
                    weight = 0.1 * (timeWt + fracWt + ((refl[j,k]/refMax)*reflWt) + ((speed[j,k]/spdMax)*windWt))
                    slip[j,k] = weight
    
    return slip

# Open the NetCDF file
inDir = "../WRF/run/"
inFile = inDir + "wrfout_d02_2022-01-16_00:00:00"
ncFile = Dataset(inFile)

# Get the variable of interest and timesteps of interest
t2 = getvar(ncFile, 'T2', timeidx=ALL_TIMES, method='cat')
dbz = getvar(ncFile, "dbz", timeidx=ALL_TIMES, method='cat')
sr = getvar(ncFile, "SR", timeidx=ALL_TIMES, method='cat')
times = list(dbz.coords['Time'].values)
umet, vmet = getvar(ncFile, "uvmet10", timeidx=ALL_TIMES, method='cat')

# Get the latitude and longitude points
lats, lons = latlon_coords(dbz)

# Get the cartopy mapping object
cart_proj = get_cartopy(dbz)
           

#### Loop through each timestep and create a new plot for each one
for i in range(len(dbz)):
    
    # The first timestep is only good for setting the initialization time for 
    # the plot title. Afterwards, only plot the subsequent timesteps.
    if i == 0:
        init = str(times[i])[11:13] + 'z ' + str(times[i])[0:10]
        
    else:
        # Create a figure
        fig = plt.figure(figsize=(20,12))
        
        # Set the GeoAxes to the projection used by WRF
        ax = plt.axes(projection=cart_proj)
    
        # Add the states and coastlines
        ax.add_feature(cfeature.STATES, edgecolor='black', linewidth=0.5)
        ax.add_feature(cfeature.COASTLINE, linewidth=0.8)
        ax.coastlines('50m', linewidth=0.8)
    
        
        ## Loop through the grid to determine the ptype at each grid location
        asi = slipIndex(times[i],t2[i],sr[i],dbz[i][0],umet[i],vmet[i])
        
        # Make filled contour for the dBZ values for rain, mix, and snow
        slipRng = np.arange(0,1.01,0.05)
        slipMesh = plt.contourf(to_np(lons), to_np(lats), to_np(asi), slipRng,
                     transform=crs.PlateCarree(), zorder=0, cmap=get_cmap("jet"))
    
    
        # Add a color bar
        cax = plt.colorbar(slipMesh, shrink=.65)
        cax.set_label('ASI')
    
        # Set the map bounds
        ax.set_xlim(cartopy_xlim(dbz))
        ax.set_ylim(cartopy_ylim(dbz))
    
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
        title_line1 = str('WRF-3km Austin Slip Index (ASI)')
        title_line2 = str('Initialized: '+init+'      Time Step: '+ step +'      Valid at: '+ hour + ' ' + date)
        plt.title(title_line1+'\n'+title_line2)
        
        # Add additional info to the plot (creator info, max/min values, etc.)
        ax.text(0,-0.01,'Matthew Toadvine\nTwitter: @toadwx', horizontalalignment='left',
                verticalalignment='top', transform=ax.transAxes)
        #ax.text(1,-0.01,'Maximum 2m Temperature: {:5.1f}\N{DEGREE SIGN}F\nMinimum Pressure: {:6.1f} hPa'.format(maxTemp,minMSLP), 
        #    horizontalalignment='right', verticalalignment='top', transform=ax.transAxes)
        
        # Save the figures
        plt.savefig('./images/asi/asi_step_'+step+'.png', format='png',bbox_inches='tight')
        #plt.show()
        plt.close()
