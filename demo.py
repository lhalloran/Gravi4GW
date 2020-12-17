import Gravi4GW
import numpy as np
import pandas as pd

# define path to DEM or GW elevation model (must be in a meter-based projection)
DEM_path = 'Geotiff/example_DEM.tif'

GW_d = 4 # assumed depth to water table from ground surface
stn_x_array = 2606457 + 20*np.arange(-10,10)
stn_y_array = 1116119 + 20*np.arange(-10,10)

file_out = 'Output/data_out.csv'

output = Gravi4GW.Gravi4GW(DEM_path, stn_x_array, stn_y_array, GW_d, accept_resid=0.02, n_r=40, do_figs=True)

dfcols = ['x [m]','y [m]','z [m]','stn. height above GW table [m]','dg/dH (x-component) [uGal/m]','dg/dH (y-component) [uGal/m]','dg/dH (z-component) [uGal/m]','dg/dH [uGal/m]']
outputdf = pd.DataFrame(data=output,columns=dfcols)
outputdf.to_csv(file_out,index=False)