import Gravi4GW
import time

DEM_path = 'geotiff/example_DEM.tif'
GW_d = 2.5 # assumed depth to water table from ground surface
stn_x_array = 2606292# + 50*np.arange(-5,5)
stn_y_array = 1116278# + 50*np.arange(-5,5)

t0=time.time()
output = Gravi4GW.Gravi4GW(DEM_path, stn_x_array, stn_y_array, GW_d, accept_resid=0.02, n_r=40, do_figs=True)
t1=time.time()
print(t1-t0)