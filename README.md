# Gravi4GW
Landon Halloran (www.ljsh.ca) 2021

## What it is and what it does.
Gravi4GW (pronounced *Gra-vee-for-ground-wa-ter*) is a python tool that enables the calculation of the conversion factor between changes in gravity (&Delta;g) as measured using time-lapse gravimetry and changes in groundwater storage (GWSC). The tool calculates &beta; (or rather, &beta;<sub>z</sub>), the rate of change in gravity as groundwater storage changes (dg/dh) in units of &mu;Gal/m<sub>H2O</sub> (=1 x 10<sup>-8</sup> s<sup>-2</sup>). A very similar quantity, wherein the effective depth to the groundwater table is set to 0, is also known as "topographical admittance.".

Basic intended uses are:
- Conversion of measured &Delta;g data to GWSC in terms of equivalent free water column height in meters.
- Creation of maps of &beta; to target gravimetric field work.
- Uncertainty analysis in hydrogeological time-lapse gravimetry studies.

![](/Output/example_output_fig.png "Example output")

## How to cite.
Halloran, L.J.S. (2022). "Improving groundwater storage change estimates using time-lapse gravimetry with Gravi4GW", *Environmental Modelling and Software* **150**, 105340. https://doi.org/10.1016/j.envsoft.2022.105340

Please contact me for more info.

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
