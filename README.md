# Gravi4GW
by Landon Halloran (www.ljsh.ca) 2021

## Overview
Gravi4GW (pronounced Gra-vee-for-ground-water) is a python tool that enables the calculation of the conversion factor between changes in gravity (&Delta;g) and changes in groundwater storage (GWSC). The tool accepts topographic or groundwater table data in the geotiff format. 

## How to cite.
The paper is:
Halloran, L.J.S. (under review) "Improving groundwater storage change estimates using time-lapse gravimetry with Gravi4GW".

## How to use. 
For a demonstration of usage, see `demo.py`.

Basic usage:
```
import Gravi4GW
# ...
# define DEM_path, stn_x_array, stn_y_array, h_eff
# ...
output = Gravi4GW.Gravi4GW(DEM_path, stn_x_array, stn_y_array, h_eff, accept_resid=0.02, n_r=40, do_figs=True)
```

