# -*- coding: utf-8 -*-
"""
Gravi4GW.py
Landon Halloran (ljsh.ca) 2021
https://github.com/lhalloran/Gravi4GW
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as interp
from osgeo import gdal

# define constants (all in SI units)
G = 6.67408E-11 # gravitational constant
rho_H2O = 1000 # density of water

def quickplotterxyz(x,y,z):
    """
    quickplotterxyz(x,y,z)
    Simple plotting of x,y,z arrays for debugging/sanity checks.
    x,y,z : 2-D arrays of x, y, and elevation
    """
    fig, axs = plt.subplots(nrows=1,ncols=3,sharex=True,sharey=True,figsize=(12,6)) 
    axs[0].imshow(x); axs[0].set_title('x, min='+str(int(min(x.flatten())))+', max='+str(int(max(x.flatten()))))
    axs[1].imshow(y); axs[1].set_title('y, min='+str(int(min(y.flatten())))+', max='+str(int(max(y.flatten()))))
    axs[2].imshow(z); axs[2].set_title('z, min='+str(int(min(z.flatten())))+', max='+str(int(max(z.flatten()))))

def pointmaker(stn_x, stn_y, max_r, nr, dens_azi=8):
    """
    pointmaker(stn_x, stn_y, max_r, nr, dens_azi=8)
    
    Returns x,y coordinates the element area represented by 
    each point. Point distribution is radial with increasing 
    distance from centre for each subsequent radius. Even 
    radial spacing for each r.
    
    Parameters:
    stn_x, stn_y : coordinates of centre point
    max_r : distance of furthest points from centre point
    nr : number of radii (including r=0) for point creation
    dens_azi : (optional) point density at first non-zero 
               radius, controls the azimuthal point 
               density. must be multiple of 4
    """
    if dens_azi%4 != 0:
        print('#Gravi4GW: pointmaker requires a radial density parameter dens_azi that is a multiple of 4')
        return 0
    rs = np.append(np.array(0.0),np.logspace(np.log10(0.1), np.log10(max_r), num=nr-1))
    log_fac = np.log10(rs[-1]/rs[-2])
    rsextra = rs[-1]*10**log_fac # for calculation of elemental areas of farthest points
    n = rs-rs
    n = dens_azi + 4*np.floor_divide(rs,10) #number of points for given r, increasing with radius
    n[0] = 1 # one point at centre
    A = rs-rs # area represent by each point at given r
    for i in np.arange(1,nr-1):
        r_out = (rs[i+1]+rs[i])/2 # outer radius of element
        r_in = (rs[i]+rs[i-1])/2 # inter radius of element
        A[i] = np.pi*(r_out**2-r_in**2)/n[i]
    A[0] = np.pi*((rs[0]+rs[1])/2)**2
    A[-1] = np.pi*((rsextra+rs[-1])**2-(rs[-1]+rs[-2])**2)/(4*n[-1])
    # make r,phi coordinates, referenced to centre point
    pt_r = np.array([])
    pt_phi = np.array([])
    pt_A = np.array([])
    for i in np.arange(nr):
        for j in np.arange(n[i]):
            pt_r = np.append(pt_r,rs[i])
            #print(i)
            #print(pt_r)
            pt_phi = np.append(pt_phi,2*np.pi*j/n[i])
            pt_A = np.append(pt_A,A[i])
    # convert to x,y cords:
    x = stn_x+np.multiply(pt_r,np.cos(pt_phi))
    y = stn_y+np.multiply(pt_r,np.sin(pt_phi))
    return x,y,pt_A

def dg(xyz_stn,xyz_pt,dm,debug=False):
    """
    dg(xyz_stn,xyz_pt,dm,debug=False) 
    
    Function that returns the delta gravity vector for a station point and 
    a grid point, given a change in mass of dm at the point
    
    xyz_stn : local coordiates (in m) of gravity station. numpy array of size 3.
    xyz_pt : local coordiates (in m) of point on GW table grid. numpy array of size 3.
    dm : (near infinitessimal) change in mass of water at the point on the grid
    debug : True will print internal values for each call of the function
    """
    dxyz = xyz_pt-xyz_stn;
    dx = dxyz[0]
    dy = dxyz[1]
    dz = dxyz[2]
    d = np.sqrt(dx**2+dy**2+dz**2) #distance between station and point
    dg1 = G*dm/d**2
    theta = np.arctan2(dz,np.sqrt(dx**2+dy**2)) # angle below xy surface
    phi = np.arctan2(dy,dx)    
    dgx = dg1*np.cos(theta)*np.sin(phi)
    dgy = dg1*np.cos(theta)*np.cos(phi)
    dgz = dg1*np.sin(theta)
    dgxyz = np.array([dgx,dgy,dgz])
    if(debug):
        print(dxyz)
        print([d,dg1,theta,phi])
        print(dgxyz)
    return dgxyz

def hillshade(array,azimuth,angle_altitude):
    # Adapted from:
    # github.com/rveciana/introduccion-python-geoespacial/blob/master/hillshade.py
    azimuth = 360.0 - azimuth 

    x, y = np.gradient(array)
    slope = np.pi/2. - np.arctan(np.sqrt(x*x + y*y))
    aspect = np.arctan2(-x, y)
    azimuthrad = azimuth*np.pi/180.
    altituderad = angle_altitude*np.pi/180.

    shaded = np.sin(altituderad)*np.sin(slope) + np.cos(altituderad)*np.cos(slope)*np.cos((azimuthrad - np.pi/2.) - aspect)

    return 255*(shaded + 1)/2


def Gravi4GW(tif_path, gravstn_x, gravstn_y, GW_depth, accept_resid=0.025, n_r=30, do_figs=True):
    """ 
    Gravi4GW.Gravi4GW()
    
    Arguments:
        tif_path: string 
            Path to the geotiff file (e.g. DEM or water table elevation model).
            Coordinate-system must be meter-based (i.e., not degrees).
        gravstn_x, gravstn_y: array-like or scalar 
            The x and y coordinate(s) of the stations in the same coordinate systems as the 
            geotiff at tif_path. If one or both of these arrays has length > 1, a grid is formed
            at all x,y pairs.
        GW_depth: scalar
            The estimated vertical distance between the gravity sensor location and the water 
            table.
     Optional arguments:
         accept_resid: approximate acceptable residual based on Bouger plate approximation (see paper)
         n_r: number of radial distances for point definition in numerical integral (default = 30)
         do_figs: boolean to enable creation of figures.
    """
    DEM_in = gdal.Open(tif_path, gdal.GA_ReadOnly) 
    print('#Gravi4GW: File '+str(tif_path)+' imported. Size = '+str(DEM_in.RasterXSize)+' x '+str(DEM_in.RasterYSize))
    DEM_z = np.array(np.float64(DEM_in.ReadAsArray()))

    # make x,y (long, lat) matrices for input DEM
    print('#Gravi4GW: Creating x,y matrices...')
    DEM_geotransform = DEM_in.GetGeoTransform()
    xy_inds = np.indices((DEM_in.RasterXSize, DEM_in.RasterYSize))
    DEM_x = DEM_geotransform[0] + DEM_geotransform[1]*xy_inds[0] + DEM_geotransform[2]*xy_inds[1]
    DEM_y = DEM_geotransform[3] + DEM_geotransform[4]*xy_inds[0] + DEM_geotransform[5]*xy_inds[1]
    DEM_x = DEM_x.transpose()
    DEM_y = DEM_y.transpose()
    if do_figs:
        quickplotterxyz(DEM_x,DEM_y,DEM_z)

    # create interpolation function
    print('#Gravi4GW: Creating DEM interpolation function...')
    interp_spline = interp.RectBivariateSpline(DEM_x[0,:],-DEM_y[:,0],DEM_z.transpose()) # x,y must be strictly increasing, hence - sign on y

    stn_array_size = [np.size(gravstn_x),np.size(gravstn_y)]
    nstns = stn_array_size[0]*stn_array_size[1]
    stn_x,stn_y = np.meshgrid(gravstn_x,gravstn_y)
    stn_x = stn_x.flatten()
    stn_y = stn_y.flatten()
    stn_z = stn_x-stn_x

    for i in np.arange(nstns):
        stn_z[i] = interp_spline(stn_x[i],-stn_y[i]) #+ 0.5
    
    max_r = GW_depth/accept_resid # max radius from stn for calculation (see notes)
    dataproc = []
    
    # cut out relevant part of DEM and x,y arrays:
    mrgn = max_r
    mmx = [np.min(gravstn_x)-mrgn, np.max(gravstn_x)+mrgn]
    mmy = [np.min(gravstn_y)-mrgn, np.max(gravstn_y)+mrgn]
    xinds = np.where(np.logical_and(np.min(DEM_x,axis=0) > mmx[0], np.max(DEM_x,axis=0) < mmx[1]))[0]
    yinds = np.where(np.logical_and(np.min(DEM_y,axis=1) > mmy[0], np.max(DEM_y,axis=1) < mmy[1]))[0]
    DEM_xC = DEM_x[yinds,:][:,xinds]
    DEM_yC = DEM_y[yinds,:][:,xinds]
    DEM_zC = DEM_z[yinds,:][:,xinds]

    for n in np.arange(nstns):
        print('Station point '+str(n+1)+' of '+str(nstns))
        stn_xyz = np.array([stn_x[n],stn_y[n],stn_z[n]])
        GW_d = GW_depth # this could be non-constant in future versions.
        
        # create sampling points and calulate DEM at these points
        print('#Gravi4GW: Defining sampling points and calculating interpolated DEM at these points...')
    
        xx,yy,AA = pointmaker(stn_xyz[0],stn_xyz[1],max_r,n_r,dens_azi=8)
        npts = np.size(xx)
        zz = xx-xx
        for i in np.arange(npts):
            zz[i] = interp_spline(xx[i],-yy[i])
        if n==0 and do_figs: # plot of points for integration (only do this once)
            plt.figure(figsize=(10,8)); plt.scatter(xx,yy,c=zz,s=AA,alpha=0.5); plt.axis('equal'); plt.title(str(npts)+' Integration Points'); plt.colorbar()# plot points
            
        # Calculate dg/dH by numerical integration:
        dH = 0.01 # drop in water table in equivalent H2O (equivalent to delta hydraulic head / porosity). Should be small so as to estimate dG/dH.
        dgxyz_sum = np.array([0,0,0])
        progresspct = 0
        print('#Gravi4GW: Evaluating delta g integral...')
        for i in np.arange(npts):
            if int(100*i/npts)-progresspct >=20: # print progress
                progresspct = int(100*i/npts)
                print('#Gravi4GW: Integration progress = '+str(progresspct)+'%')
            dm = dH*rho_H2O*AA[i]
            pt_xyz = np.array([xx[i],yy[i],zz[i]-GW_d])
            dgxyz_sum = dgxyz_sum + dg(stn_xyz,pt_xyz,dm)
        dg1 = np.append(dgxyz_sum,np.sqrt(dgxyz_sum[0]**2+dgxyz_sum[1]**2+dgxyz_sum[2]**2))
        dgdH = dg1/dH
        dgdH_uGal = dgdH*1E8 #in microGal/mH2O
        dataproc.append([stn_xyz[0],stn_xyz[1],stn_xyz[2],GW_d,dgdH_uGal[0],dgdH_uGal[1],dgdH_uGal[2],dgdH_uGal[3]])
        print('#Gravi4GW: beta_x = ' + str(dgdH_uGal[0]) + ' uGal/mH2O')
        print('#Gravi4GW: beta_y = ' + str(dgdH_uGal[1]) + ' uGal/mH2O')
        print('#Gravi4GW: beta_z = ' + str(dgdH_uGal[2]) + ' uGal/mH2O')
        print('#Gravi4GW: beta = '   + str(dgdH_uGal[3]) + ' uGal/mH2O')
    dataproc = np.array(dataproc) #convert data to numpy array

    #%% plot the results
    only1xy = np.size(gravstn_x) > 1 and np.size(gravstn_y) > 1
    if do_figs and only1xy:
        fig, axs = plt.subplots(nrows=1,ncols=2,sharex=True,sharey=True,figsize=(15,8))
        DEM_hs = hillshade(DEM_zC,45,20)
        cbobj1 = axs[0].imshow(np.flipud(DEM_hs),cmap='Greys',alpha=1, interpolation='bilinear',extent=[DEM_xC[0,0],DEM_xC[0,-1],DEM_yC[0,0],DEM_yC[-1,0]])
        demobj = axs[0].contourf(DEM_xC,DEM_yC,DEM_zC,cmap='gist_earth',alpha=0.5, levels=80)
        plt.gca().invert_yaxis()
        
        axs[0].set_title('Input data (m)')
        axs[0].set_aspect('equal')
        axs[0].grid(c='k')
        plt.setp(axs[0].get_xticklabels(), rotation=45)
        plt.setp(axs[0].get_yticklabels(), rotation=45)
        
        cbarax1 = fig.add_axes([0.48, 0.2, 0.01, 0.6])
        fig.colorbar(demobj, cax=cbarax1, orientation='vertical')
        
        levels = np.linspace(np.min(dataproc[:,-1]),41.93*2-np.min(dataproc[:,-1]),101)
        cbobj2 = axs[1].contourf(gravstn_x, gravstn_y, dataproc[:,-1].reshape(np.flip(stn_array_size)),levels=levels,cmap='bwr')
        for c in cbobj2.collections:
            c.set_edgecolor("face")
        axs[1].set_title(r'$\beta$ ($\mu$Gal/mH$_2$O)')
        axs[1].set_aspect('equal')
        axs[1].grid(c='k')
        plt.gca().invert_yaxis()
        plt.setp(axs[1].get_xticklabels(), rotation=45)
        cbarax2 = fig.add_axes([0.91, 0.2, 0.01, 0.6])
        fig.colorbar(cbobj2, cax=cbarax2, orientation='vertical')

    return dataproc
