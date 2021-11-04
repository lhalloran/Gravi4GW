# Gravi4GW
Landon Halloran (www.ljsh.ca) 2021

## What it is and what it does.
Gravi4GW (pronounced *Gra-vee-for-ground-wa-ter*) is a python tool that enables the calculation of the conversion factor between changes in gravity (&Delta;g) as measured using time-lapse gravimetry and changes in groundwater storage (GWSC). The tool calculates &beta;, the rate of change in gravity as groundwater storage changes (dg/dh) in units of &mu;Gal/m<sub>H2O</sub> (=1 x 10<sup>-8</sup> s<sup>-2</sup>). This term is also known as "topographical admittance.".

Basic intended uses are:
- Conversion of measured &Delta;g data to GWSC in terms of equivalent free water column height in meters.
- Creation of maps of &beta; to target gravimetric field work.
- Uncertainty analysis in hydrogeological time-lapse gravimetry studies.

![](/Output/example_output_fig.png "Example output")

## How to cite.
Please contact me for more info.
Halloran, L.J.S. (in revision), "Improving groundwater storage change estimates using time-lapse gravimetry with Gravi4GW".
A citable pre-print version is available here: https://doi.org/10.1002/essoar.10507438.1

## How to use. 
For a simple demonstration of usage, see `demo.py`.

Basic usage:
```
import Gravi4GW
# ...
# User defines GEOTIFF_path, stn_x, stn_y, h_eff ...
# ...
output = Gravi4GW.Gravi4GW(GEOTIFF_path, stn_x, stn_y, h_eff)
```

More information on the inputs and outputs of the functions contained in `Gravi4GW.py` can be obtained by typing, for example, `help(Gravi4GW.Gravi4GW)`.