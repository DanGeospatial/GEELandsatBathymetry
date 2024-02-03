# GEELandsatBathymetry
Script used to calculate shallow water bathymetry with Landsat imagery in Google Earth Engine.

Automatically switched between Landsat satellite based on date range. 
Requires batch export function in GEETools library developed by Rodrigo E. Principe.

Uses bottom albedo-independent Bathymetry algorithm modified by Li et al. (2021) and originally developed by Stumpf and Holderied (2003).

WARNING

The results of this tool have not been validated and may be incorrect!
You must understand how to produce a correct bathymetry estimation equation for this to produce meaningful results.

This function does not use any coefficients to standardize the images between satellite sensors.
