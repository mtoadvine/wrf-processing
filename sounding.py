#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Create Skew-T Log-P sounding for a single point of interest
from WRF output.

@author: Matthew Toadvine
@date:   01/14/2022
"""

# Import required libraries
import sys, os
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from metpy.plots import SkewT,Hodograph
import metpy.calc as mpcalc
from metpy.units import units
from PIL import Image
from wrf import (getvar,ALL_TIMES,ll_to_xy,xy_to_ll)

# Open the NetCDF file
inDir = "../WRF/run/"
inFile = inDir + "wrfout_d01_2022-03-23_12:00:00"
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

# If the Lat/Lon point was valid and converted to nearest x/y coordinates,
# convert the nearest x/y point back to Lat/Lon for title
lat, lon = xy_to_ll(ncfile,x,y)

# Get the variable of interest and timesteps of interest
# [times,vert levels,lat,lon]
press = getvar(ncfile, 'pressure', timeidx=ALL_TIMES, method='cat')[:,:,y,x]
height = getvar(ncfile, 'height_agl', timeidx=ALL_TIMES, method='cat')[:,:,y,x]
temp = getvar(ncfile, 'temp', timeidx=ALL_TIMES, method='cat')[:,:,y,x]
dewp = getvar(ncfile, "td", timeidx=ALL_TIMES, method='cat')[:,:,y,x]
umet, vmet = getvar(ncfile, "uvmet", timeidx=ALL_TIMES, method='cat')
umet = umet[:,:,y,x]
vmet = vmet[:,:,y,x]
times = list(temp.coords['Time'].values)

# Define a list to contain all of the image files for automatic 
# GIF processing and remove all old images from the directory
frames = []
os.system("rm ./images/sounding/*")

#### Loop through each timestep and create a new plot for each one
for i in range(len(times)):
    
    # The first timestep is only good for setting the initialization time for 
    # the plot title. Afterwards, only plot the subsequent timesteps.
    if i == 0:
        init = str(times[i])[11:13] + 'z ' + str(times[i])[0:10]
        
    else:
        
        # Create numpy arrays to hold variable values for metpy use
        p = np.array(press[i].values)*units.hPa
        z = np.array(height[i].values)*units.m
        t = np.array(temp[i].values)*units('K')
        td = np.array(dewp[i].values)*units.degC
        u = np.array(umet[i].values)*units('m/s')
        v = np.array(vmet[i].values)*units('m/s')
        
        # Convert temp to deg C from Kelvin
        t = t.to('degC')
        
        # Convert u and v components to knots
        u = u.to('knots')
        v = v.to('knots')
        
        ##### Print Skew-T Plots (a lot used from Sarah Purpura) #####

        # Change the DPI of the resulting figure. Higher DPI drastically improves the
        # look of the text rendering.
        plt.rcParams['savefig.dpi'] = 300
        
        #print()
        ## Change default to be better for skew-T
        fig = plt.figure(figsize=(12,10))
        skew = SkewT(fig,rotation=45)

        ## Plot the data using normal plotting functions, in this case using
        ## log scaling in Y, as dictated by the typical meteorological plot
        skew.plot(p, t, 'r')
        skew.plot(p, td, 'g')
        skew.plot_barbs(p[:-4:2], u[:-4:2], v[:-4:2])

        ## Add the relevant special lines
        skew.plot_dry_adiabats(t0=np.arange(233,533,10)*units('K'),alpha=0.25)
        skew.plot_moist_adiabats(t0=np.arange(233,323,5)*units('K'),alpha=0.25)
        skew.plot_mixing_lines(alpha=0.25)
        skew.ax.set_ylim(1050, 100)
        skew.ax.set_xlim(-30,40)

        ## Add parcel profile to plot as black line
        ## Requires that the variables have associated units
        prof = mpcalc.parcel_profile(p, t[0], td[0])
        skew.plot(p, prof, 'k', linewidth=2,alpha=0.25)
        #skew.shade_cin(p,t,prof,td)
        #skew.shade_cape(p,t,prof)
        
        '''
        ~~~~~~ I have not adjusted this to match the current data structures
        ~~~~~~ so variables may need to be corrected before use
        ## This will add a line and place where the LCL is located
        #skew.plot(lcl_p,lcl_t.to('degC'), marker='_',  markersize=22, color='blue')
        #plt.text(lcl_t.to('degC').m+5,lcl_p.m+20,'LCL',color='blue',size=14)
        
        ##this will add a line and place where the LFC is located
        skew.plot(lfc_p,lfc_t.to('degC'), marker='_',  markersize=22, color='blue')
        plt.text(lfc_t.to('degC').m+5,lfc_p.m+20,'LFC',color='blue',size=14)

        ##this will add a line and place where the EL is located
        skew.plot(el_p,el_t.to('degC'), marker='_',  markersize=22, color='blue')
        plt.text(el_t.to('degC').m+5,el_p.m+20,'EL',color='blue',size=14)
        '''
        
        
        #### Add plot titles and basic model run information
        
        # Define timestep hour and date for title    
        date = str(times[i])[0:10]
        hour = str(times[i])[11:16] + 'z'
        if i < 10:
            step = '0' + str(i)
        else:
            step = str(i)
        
        # Create two title lines
        title_line1 = str('WRF-3km Sounding for ' + str(np.around(lat.values,2)) + 'N, ' + str(np.around(-1*lon.values,2)) + 'W')
        title_line2 = str('Initialized: '+init+'      Time Step: '+ step +'      Valid at: '+ hour + ' ' + date)
        plt.title(title_line1+'\n'+title_line2)
        
        ## Make my own hodograph from scratch
        # Weird placement, but placing this block above the title causes 
        # the title to align based off of the hodograph, not the skewt
        ax_hod = inset_axes(skew.ax, '25%', '25%', loc=1)
        ho = Hodograph(ax_hod, component_range=80.)
        ho.add_grid(increment=20)
        ho.plot_colormapped(u[p >= 100.0 * units.hPa], 
                            v[p >= 100.0 * units.hPa], 
                            z[p >= 100.0 * units.hPa])  # Plot a line colored by wind speed
        
        # Add additional info to the plot (creator info, max/min values, etc.)
        skew.ax.text(0,-0.05,'Matthew Toadvine\nTwitter: @toadwx', horizontalalignment='left',
                verticalalignment='top', transform=skew.ax.transAxes)

        
        ## Show the plot
        plt.savefig('images/sounding/sounding'+step+'.png',bbox_inches='tight')
        #plt.show()
        plt.close()
        
        # Append the image to the list frames for GIF processing
        image = './images/sounding/sounding'+step+'.png'
        frames.append(Image.open(image))

# Perform automatic GIF processing and save the file
frame_one = frames[0]
frame_one.save("./images/sounding/sounding.gif", format="GIF", append_images=frames,
           save_all=True, duration=500, loop=0)
