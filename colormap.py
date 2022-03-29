#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Create color maps for WRF output plotting.

@author: Matthew Toadvine
@date:   01/14/2022
"""

# Import matplotlib library
import matplotlib

# National Weather Service QPF
nws_precip_colors = [
    "#04e9e7",  # 0.01 - 0.10 inches
    "#019ff4",  # 0.10 - 0.25 inches
    "#0300f4",  # 0.25 - 0.50 inches
    "#02fd02",  # 0.50 - 0.75 inches
    "#01c501",  # 0.75 - 1.00 inches
    "#008e00",  # 1.00 - 1.50 inches
    "#fdf802",  # 1.50 - 2.00 inches
    "#e5bc00",  # 2.00 - 2.50 inches
    "#fd9500",  # 2.50 - 3.00 inches
    "#fd0000",  # 3.00 - 4.00 inches
    "#d40000",  # 4.00 - 5.00 inches
    "#bc0000",  # 5.00 - 6.00 inches
    "#f800fd",  # 6.00 - 8.00 inches
    "#9854c6",  # 8.00 - 10.00 inches
    "#fdfdfd"]   # 10.00+

precip_colormap = matplotlib.colors.ListedColormap(nws_precip_colors)

# National Weather Service Reflectivity
nws_dbz_colors = [
    "#FFFFFF", # 0
    "#808080", # 5
    "#00FFFF", # 5-10
    "#00BFFF", # 10-15
    "#0000FF", # 15-20
    "#00FF00", # 20-25
    "#32CD32", # 25-30
    "#008000", # 30-35
    "#FFFF00", # 35-40
    "#DAA520", # 40-45
    "#FFA500", # 45-50
    "#FF0000", # 50-55
    "#8B0000", # 55-60
    "#800000", # 60-65
    "#FF00FF", # 65-70
    "#8A2BE2", # 70-75
    "#FFFFFF"] # 75+

dbz_colormap = matplotlib.colors.ListedColormap(nws_dbz_colors)

# Snow Reflectivity
snow_cmap = [
    "#bfe3f6",
    "#aed4eb",
    "#9ec5df",
    "#8eb7d4",
    "#7ea8c9",
    "#709abe",
    "#628bb3",
    "#547da8",
    "#486f9d",
    "#3c6192",
    "#305487",
    "#25467b",
    "#1a3970",
    "#0f2d64",
    "#042058"]

snow_colormap = matplotlib.colors.ListedColormap(snow_cmap)


# Sleet Reflectivity
sleet_cmap = [
    "#f6bfe9",
    "#eab2dd",
    "#dea4d1",
    "#d397c5",
    "#c78ab9",
    "#bc7dad",
    "#b070a2",
    "#a56496",
    "#9a578b",
    "#8f4a80",
    "#843e75",
    "#79316a",
    "#6e2460",
    "#631655",
    "#58044b"]

sleet_colormap = matplotlib.colors.ListedColormap(sleet_cmap)

# Freezing Rain Reflectivity
zr_cmap = [
    "#f7d6be",
	"#f5ceb1",
	"#f3c5a4",
	"#f1bd97",
	"#efb48a",
	"#ecac7d",
	"#e9a371",
	"#e79a65",
	"#e49259",
	"#e0894d",
	"#dd8041",
	"#d97734",
	"#d66f28",
	"#d26519",
	"#ce5c06"]

zr_colormap = matplotlib.colors.ListedColormap(zr_cmap)

# Snow Accumulation
snowAcc_cmap = [
    "#a8d9eb",  # 0-1"
    "#91d1f4",  # 1-2"
    "#82c8fd",  # 2-3"
    "#7ebcff",  # 3-4"
    "#89aeff",  # 4-5"
    "#9f9dff",  # 5-6"
    "#ba88fa",  # 6-7"
    "#d66ee7"]  # 7-8"

snowAcc_colormap = matplotlib.colors.ListedColormap(snowAcc_cmap)