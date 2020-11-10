import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as interp
from osgeo import gdal
import Gravi4GW, G4GW_f

# import DEM (must be in a meter-based projection)
#DEM_path = 'DEMs/Tsalet_30m_DEM_EPSG2056.tif'
DEM_path = 'DEMs/RechyDEM_swisstopo_2m.tif'
#DEM_path = 'DEMs/Salar_de_Uyuni_DEM.tif'
#DEM_path = 'DEMs/Maules_Creek_20m_DEM-H.tif'
# define x,y of gravity station
GW_d = 15 # assumed depth to water table from ground surface
stn_x_array = 2606292 + 50*np.arange(-5,5)
stn_y_array = 1116278 + 50*np.arange(-5,5)
#stn_x_array = 1327850 + 30*np.arange(-20,20)
#stn_y_array = -2236035.0 + 30*np.arange(-20,20)
#stn_x_array = 1508000+ 20*np.arange(-20,20)
#stn_y_array = -3523500.0+ 20*np.arange(-20,20)

output = Gravi4GW.Gravi4GW(DEM_path, stn_x_array, stn_y_array, GW_d, accept_resid=0.025, n_r=40, single_or_array='array')
