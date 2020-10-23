# -*- coding: utf-8 -*-
"""
Gravi4GW###.py
14.10.2020 L Halloran

Part of time-lapse gravimetry project. Determines change in gravity caused by 
drop/rise of the water table, taking topography into account.

v001/002: created. stumbling upon efficiently making fine grid
v003/004: fine grid idea abandoned, now making concentric points array, and interpolating that way
v005: moved functions to another file G4GW_f... might convert to class later
v006: now synced to github repo Gravi4GW
"""
import numpy as np
#import pandas as pd
import matplotlib.pyplot as plt
import scipy.interpolate as interp
from osgeo import gdal

import G4GW_f # could turn this into a class-style python package...

# import DEM (must be in a meter-based projection)
#DEM_path = 'DEMs/Tsalet_30m_DEM_EPSG2056.tif'
DEM_path = 'DEMs/RechyDEM_swisstopo_2m.tif'
#DEM_path = 'DEMs/Salar_de_Uyuni_DEM.tif'
DEM_in = gdal.Open(DEM_path, gdal.GA_ReadOnly) 
print('#Gravi4GW: DEM file '+str(DEM_path)+' imported. Size = '+str(DEM_in.RasterXSize)+' x '+str(DEM_in.RasterYSize))
DEM_z = np.array(np.float64(DEM_in.ReadAsArray()))

#for debugging:
#DEM_z = DEM_z-DEM_z+ 2437.7

# make x,y (long, lat) matrices for input DEM
print('#Gravi4GW: Creating x,y matrices...')
DEM_geotransform = DEM_in.GetGeoTransform()
xy_inds = np.indices((DEM_in.RasterXSize, DEM_in.RasterYSize))
DEM_x = DEM_geotransform[0] + DEM_geotransform[1]*xy_inds[0] + DEM_geotransform[2]*xy_inds[1]
DEM_y = DEM_geotransform[3] + DEM_geotransform[4]*xy_inds[0] + DEM_geotransform[5]*xy_inds[1]
DEM_x=DEM_x.transpose()
DEM_y=DEM_y.transpose()
G4GW_f.quickplotterxyz(DEM_x,DEM_y,DEM_z)

# create interpolation function
print('#Gravi4GW: Creating DEM interpolation function...')
# ...TO DO HERE???: cut out portion of input data that is relevant to points (decrease size of input to RectBivariateSpline)
interp_spline = interp.RectBivariateSpline(DEM_x[0,:],-DEM_y[:,0],DEM_z.transpose()) # x,y must be strictly increasing, hence - sign on y

# define x,y of gravity station
GW_d = 4 # assumed depth to water table from ground surface
stn_x_array = 2606190.1 + 20*np.arange(-10,10)
stn_y_array = 1116135.25 + 20*np.arange(-10,10)
stn_array_size = [np.size(stn_x_array),np.size(stn_y_array)]
nstns = stn_array_size[0]*stn_array_size[1]
stn_x,stn_y = np.meshgrid(stn_x_array,stn_y_array)
stn_x = stn_x.flatten()
stn_y = stn_y.flatten()
stn_z = stn_x-stn_x

# cut out relevant part of DEM and x,y arrays:
mrgn = 200
mmx = [np.min(stn_x_array)-mrgn, np.max(stn_x_array)+mrgn]
mmy = [np.min(stn_y_array)-mrgn, np.max(stn_y_array)+mrgn]
xinds = np.where(np.logical_and(np.min(DEM_x,axis=0) > mmx[0], np.max(DEM_x,axis=0) < mmx[1]))[0]
yinds = np.where(np.logical_and(np.min(DEM_y,axis=1) > mmy[0], np.max(DEM_y,axis=1) < mmy[1]))[0]
DEM_xC = DEM_x[yinds,:][:,xinds]
DEM_yC = DEM_y[yinds,:][:,xinds]
DEM_zC = DEM_z[yinds,:][:,xinds]

# for 
for i in np.arange(nstns):
    stn_z[i] = interp_spline(stn_x[i],-stn_y[i]) + 0.5

dataproc = []
for n in np.arange(nstns):
    print(str(n))
    stn_xyz=np.array([stn_x[n],stn_y[n],stn_z[n]])
    GW_d = GW_d # this could be non-constant in future.
    
    # create sampling points and calulate DEM at these points
    print('#Gravi4GW: Defining sampling points and calculating interpolated DEM at these points...')
    n_r = 40 # number of radial distances for point definition
    resid_OK = 0.05 # approximate acceptable residual based on Bouger plate approx. (see notes)
    max_r = GW_d/resid_OK # max radius from stn for calculation (see notes)
    xx,yy,AA= G4GW_f.point_maker(stn_xyz[0],stn_xyz[1],max_r,n_r,dens_rad=8)
    npts=np.size(xx)
    zz=xx-xx
    for i in np.arange(npts):
        zz[i]=interp_spline(xx[i],-yy[i])
    if n==0: # plot of points for integration (only do this once)
        plt.figure(figsize=(10,8)); plt.scatter(xx,yy,c=zz,s=AA,alpha=0.5); plt.axis('equal'); plt.title(str(npts)+' Integration Points'); plt.colorbar()# plot points
    
    # Calculate dg/dH by numerical integration:
    dH=0.01 # drop in water table in equivalent H2O (equivalent to delta hydraulic head / porsity). Should be small so as to estimate dG/dH.
    dgxyz_sum = np.array([0,0,0])
    progresspct=0
    print('#Gravi4GW: Evaluating delta g integral...')
    for i in np.arange(npts):
        if int(100*i/npts)-progresspct >=10: # print progress
            progresspct=int(100*i/npts)
            print('#Gravi4GW: Integration progress = '+str(progresspct)+'%')
        dm = dH*G4GW_f.rho_H2O*AA[i]
        pt_xyz = np.array([xx[i],yy[i],zz[i]-GW_d])
        dgxyz_sum = dgxyz_sum + G4GW_f.dg(stn_xyz,pt_xyz,dm)
    dg = np.append(dgxyz_sum,np.sqrt(dgxyz_sum[0]**2+dgxyz_sum[1]**2+dgxyz_sum[2]**2))
    dgdH = dg/dH
    dgdH_uGal = dgdH*1E8 #in microGal/mH2O
    dataproc.append([stn_xyz[0],stn_xyz[1],stn_xyz[2],GW_d,dgdH_uGal[0],dgdH_uGal[1],dgdH_uGal[2],dgdH_uGal[3]])
    print('#Gravi4GW: dg_x = ' + str(dgdH_uGal[0]) + ' uGal/mH2O')
    print('#Gravi4GW: dg_y = ' + str(dgdH_uGal[1]) + ' uGal/mH2O')
    print('#Gravi4GW: dg_z = ' + str(dgdH_uGal[2]) + ' uGal/mH2O')
    print('#Gravi4GW: dg = '   + str(dgdH_uGal[3]) + ' uGal/mH2O')
dataproc=np.array(dataproc) #convert data to numpy array

# maybe improve this using https://matplotlib.org/3.1.1/gallery/specialty_plots/topographic_hillshading.html
fig, axs = plt.subplots(nrows=1,ncols=2,sharex=True,sharey=True,figsize=(12,6)) 
DEM_hs = G4GW_f.hillshade(DEM_zC,45,30)
demobj = axs[0].contourf(DEM_xC,DEM_yC,DEM_zC,cmap='gist_earth',alpha=1, levels=15)
axs[0].contourf(DEM_xC,DEM_yC,DEM_hs,cmap='Greys',alpha=0.65,levels=20)
axs[0].set_title('Input DEM')
cbarax1 = fig.add_axes([0.48, 0.2, 0.01, 0.6])
fig.colorbar(demobj, cax=cbarax1, orientation='vertical')

cbobj = axs[1].contourf(stn_x_array, stn_y_array, dataproc[:,-1].reshape(stn_array_size))
axs[1].set_title('dg/dH (uGal/mH2O)')
cbarax2 = fig.add_axes([0.91, 0.2, 0.01, 0.6])
fig.colorbar(cbobj, cax=cbarax2, orientation='vertical')
