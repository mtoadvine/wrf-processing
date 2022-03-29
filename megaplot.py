#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The MEGA-PLOTTER SCRIPT

The aim of this script is to generate several types of plots and
(hopefully) do so quickly with mulitprocessing. The types of plots
generated by this script are listed below.

Plots:
    - 2m Temperature
    - QPF
    - Sounding at a single lat/lon
    - Precipitation Type Reflectivity
    - Standard Reflectivity
    - Snow/Ice Accumulation
    - Austin Slip Index (ASI) [a joke worth plotting]

@author:        Matthew Toadvine
@version:       1.0
@edited:        02/04/2022
"""

# Import required libraries
import sys
import multiprocessing as mp
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import colormap as cm
import cartopy.crs as crs
import cartopy.feature as cfeature
import metpy.calc as mpcalc
import metpy.units as units
from math import *
from wrf import (to_np, getvar, smooth2d, get_cartopy, cartopy_xlim,
                 cartopy_ylim, latlon_coords, ALL_TIMES, ll_to_xy)


##############################################################################
# ASSIGN GLOBAL VARIABLES                                                    #
#                                                                            #
# These are some basic assignments that will be used throughout the entirety #
# of the script                                                              #
##############################################################################

# Open the NetCDF file
inDir = "../WRF/run/"
inFile = inDir + "wrfout_d01_2022-02-07_00:00:00"
ncfile = Dataset(inFile)

# Set lat/lon for point data (sounding, temps, etc)
lat = 35.1919
lon = -80.6827

# Check to see if the Lat/Lon point is within the domain
# If invalid, quit with error message
try:
    x, y = ll_to_xy(ncfile,lat,lon)
except:
    sys.exit('Lat/Lon point is not within the domain.')


##############################################################################
#   2 METER TEMPERATURE                                                      #
##############################################################################

def t2m(inFile,longitude,latitude):
    
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
    
    
    #### Loop through each timestep and create a new plot for each timestep
    for i in range(len(temp2m)):
        
        # The first timestep is only good for setting the initialization time for 
        # the plot title. Afterwards, only plot the subsequent timesteps.
        if i == 0:
            init = str(times[i])[11:13] + 'z ' + str(times[i])[0:10]
            
        else:
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
                         transform=crs.PlateCarree(), zorder=0, cmap=get_cmap("jet"))
            frzLine = plt.contour(to_np(lons), to_np(lats), to_np(temp2m[i]), np.arange(32,33),
                                   transform=crs.PlateCarree(), zorder=0, colors="black", linestyles='dashed')
        
            # Add a color bar
            cax = plt.colorbar(tempMesh, shrink=.98)
            cax.set_label('\N{DEGREE SIGN}F')
        
            # Set the map boundaries based on data extent
            ax.set_xlim(cartopy_xlim(smooth_slp))
            ax.set_ylim(cartopy_ylim(smooth_slp))
        
            # Add the gridlines
            ax.gridlines(color="black", linestyle="dotted")
            
            # Find minimum MSLP value for display
            #minTemp = np.amin(temp2m[i].squeeze())
            #minTemp = minTemp.values
        
            # Find maximum temperature for display
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


##############################################################################
#   QPF                                                                      #
##############################################################################

def qpf(inFile):
    


##############################################################################
#   SKEW-T LOG-P                                                             #
##############################################################################

def skewT(inFile,longitude,latitude):



##############################################################################
#   PRECIPITATION-TYPE REFLECTIVITY                                          #
##############################################################################

def ptype(inFile):
    


##############################################################################
#   REFLECTIVITY (STANDARD)                                                  #
##############################################################################

def refl(inFile):
    
   
    
##############################################################################
#   SNOW AND ICE ACCUMULATION                                                #
##############################################################################

def snowAccum(inFile):

    

##############################################################################
#   AUSTIN SLIP INDEX (ASI)                                                  #
##############################################################################

def asi(inFile):
    
    
    
    
    
##############################################################################
#   MAIN PROGRAM                                                             #
#                                                                            #
#   Here is where the main program is run. Multiprocessing is (hopefully)    #
#   going to be used as the above functions are called to speed up the       #
#   processing time.                                                         #
##############################################################################

if __init__() == __main__():
    