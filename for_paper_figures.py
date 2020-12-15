import Gravi4GW
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from osgeo import gdal
from matplotlib import cm
#import matplotlib.patches as mpatches

#%%
DEM_path = 'DEMs/RechyDEM_swisstopo_2m.tif'
stn_x_array = 2606457 + 20*np.arange(-25,25)
stn_y_array = 1116119 + 20*np.arange(-25,25)
dfcols = ['x [m]','y [m]','z [m]','stn. height above GW table [m]','dg/dH (x-component) [uGal/m]','dg/dH (y-component) [uGal/m]','dg/dH (z-component) [uGal/m]','dg/dH [uGal/m]']
#%% decomment this to run first time (if no csv output files exist)
#GW_d = 2.5 # assumed depth to water table from ground surface
#output1 = Gravi4GW.Gravi4GW(DEM_path, stn_x_array, stn_y_array, GW_d, accept_resid=0.02, n_r=40, do_figs=True)
file_out1 = 'Output/data_out250cm.csv'
#outputdf1 = pd.DataFrame(data=output1,columns=dfcols)
#outputdf1.to_csv(file_out1,index=False)
#
#GW_d = 5 # assumed depth to water table from ground surface
#output2 = Gravi4GW.Gravi4GW(DEM_path, stn_x_array, stn_y_array, GW_d, accept_resid=0.02, n_r=40, do_figs=True)
file_out2 = 'Output/data_out500cm.csv'
#outputdf2 = pd.DataFrame(data=output2,columns=dfcols)
#outputdf2.to_csv(file_out,index=False)
#
#GW_d = 7.5 # assumed depth to water table from ground surface
#output3 = Gravi4GW.Gravi4GW(DEM_path, stn_x_array, stn_y_array, GW_d, accept_resid=0.02, n_r=40, do_figs=True)
file_out3 = 'Output/data_out750cm.csv'
#outputdf3 = pd.DataFrame(data=output3,columns=dfcols)
#outputdf3.to_csv(file_out,index=False)

output1=np.array(pd.read_csv(file_out1))
output2=np.array(pd.read_csv(file_out2))
output3=np.array(pd.read_csv(file_out3))

#%%
DEM_in = gdal.Open(DEM_path, gdal.GA_ReadOnly) 
DEM_z = np.array(np.float64(DEM_in.ReadAsArray()))
DEM_geotransform = DEM_in.GetGeoTransform()
xy_inds = np.indices((DEM_in.RasterXSize, DEM_in.RasterYSize))
DEM_x = DEM_geotransform[0] + DEM_geotransform[1]*xy_inds[0] + DEM_geotransform[2]*xy_inds[1]
DEM_y = DEM_geotransform[3] + DEM_geotransform[4]*xy_inds[0] + DEM_geotransform[5]*xy_inds[1]
DEM_x=DEM_x.transpose()
DEM_y=DEM_y.transpose()
mmx=[np.min(output1[:,0]),np.max(output1[:,0])]
mmy=[np.min(output1[:,1]),np.max(output1[:,1])]
xinds = np.where(np.logical_and(np.min(DEM_x,axis=0) > mmx[0], np.max(DEM_x,axis=0) < mmx[1]))[0]
yinds = np.where(np.logical_and(np.min(DEM_y,axis=1) > mmy[0], np.max(DEM_y,axis=1) < mmy[1]))[0]
DEM_xC = DEM_x[yinds,:][:,xinds]
DEM_yC = DEM_y[yinds,:][:,xinds]
DEM_zC = DEM_z[yinds,:][:,xinds]

#%%
fig, axs = plt.subplots(nrows=4,ncols=1,sharex=True,sharey=True,figsize=(6,16))

DEM_hs = Gravi4GW.hillshade(DEM_zC,45,20)
#ls = LightSource(azdeg=315, altdeg=45)
cbobj1 = axs[0].imshow(np.flipud(DEM_hs),cmap='Greys',alpha=1, interpolation='bilinear',extent=[DEM_xC[0,0],DEM_xC[0,-1],DEM_yC[0,0],DEM_yC[-1,0]])
demobj = axs[0].contourf(DEM_xC,DEM_yC,DEM_zC,cmap='gist_earth',alpha=0.5, levels=80)
plt.gca().invert_yaxis()

axs[0].set_title('a)')
axs[0].set_aspect('equal')
axs[0].grid(c='k')
plt.setp(axs[0].get_xticklabels(), rotation=45)
plt.setp(axs[0].get_yticklabels(), rotation=45)

cbarax1 = fig.add_axes([0.80, 0.75, 0.04, 0.2])
fig.colorbar(demobj, cax=cbarax1, orientation='vertical')
levels=np.linspace(np.min(output1[:,-1]),41.93*2-np.min(output1[:,-1]),101)
#cbobj2 = axs[1].contourf(stn_x_array, stn_y_array, dataproc[:,-1].reshape(stn_array_size),levels=levels,cmap='bwr')
stn_array_size = [np.size(stn_x_array),np.size(stn_y_array)]
cbobj1 = axs[1].contourf(output1[:,0].reshape(np.flip(stn_array_size)), output1[:,1].reshape(np.flip(stn_array_size)), output1[:,-1].reshape(np.flip(stn_array_size)),levels=levels,cmap='bwr')
for c in cbobj1.collections:
    c.set_edgecolor("face")
axs[1].set_title('b)')
axs[1].set_aspect('equal')
axs[1].grid(c='k')
#plt.gca().invert_yaxis()
plt.setp(axs[1].get_yticklabels(), rotation=45)
plt.setp(axs[1].get_xticklabels(), rotation=45)
#cbarax1 = fig.add_axes([0.91, 0.2, 0.01, 0.6])
#fig.colorbar(cbobj1, cax=cbarax1, orientation='vertical')

cbobj2 = axs[2].contourf(output2[:,0].reshape(np.flip(stn_array_size)), output2[:,1].reshape(np.flip(stn_array_size)), output2[:,-1].reshape(np.flip(stn_array_size)),levels=levels,cmap='bwr')
for c in cbobj2.collections:
    c.set_edgecolor("face")
axs[2].set_title('c)')
axs[2].set_aspect('equal')
axs[2].grid(c='k')
#plt.gca().invert_yaxis()
plt.setp(axs[2].get_yticklabels(), rotation=45)
plt.setp(axs[2].get_xticklabels(), rotation=45)
#cbarax2 = fig.add_axes([0.91, 0.2, 0.01, 0.6])
#fig.colorbar(cbobj2, cax=cbarax2, orientation='vertical')

cbobj3 = axs[3].contourf(output3[:,0].reshape(np.flip(stn_array_size)), output3[:,1].reshape(np.flip(stn_array_size)), output3[:,-1].reshape(np.flip(stn_array_size)),levels=levels,cmap='bwr')
for c in cbobj3.collections:
    c.set_edgecolor("face")
axs[3].set_title('d)')
axs[3].set_aspect('equal')
axs[3].grid(c='k')
#plt.gca().invert_yaxis()
plt.setp(axs[3].get_yticklabels(), rotation=45)
plt.setp(axs[3].get_xticklabels(), rotation=45)
cbarax3 = fig.add_axes([0.80, 0.4, 0.04, 0.2])
fig.colorbar(cbobj3, cax=cbarax3, orientation='vertical')


#plt.savefig('Output/Rechy_3depths3.pdf')
#%% influence of depth...
n_stns = 21
stn_xa = np.linspace(2606277,2606657,n_stns)
stn_ya = np.linspace(1116299,1115959,n_stns)

# calculate beta at all depths, for all defined stns (slightly inefficient...)
incdpthout=[]
GW_da = np.arange(2,32,2)
for i in np.arange(n_stns):
    print('stn = '+str(i))
    incdpthout1=[]
    stn_x=stn_xa[i]
    stn_y=stn_ya[i]
    for d in GW_da:
        out = Gravi4GW.Gravi4GW(DEM_path, stn_x, stn_y, d, accept_resid=0.02, n_r=40, do_figs=False)
        incdpthout1.append(out)
    incdpthout1 = np.array(incdpthout1).squeeze()
    incdpthout.append(incdpthout1)
incdpthout = np.array(incdpthout).squeeze()

#%% set up for plots
colours = cm.get_cmap('cividis_r',n_stns)
fig, axs = plt.subplots(nrows=2,ncols=2,sharex=False,sharey=False,figsize=(10,10))

# make map with selected points overlayed in corresponding colours
#plt.figure()
#fig,ax = plt.subplots()
DEM_hs = Gravi4GW.hillshade(DEM_zC,45,20)
#ls = LightSource(azdeg=315, altdeg=45)
cbobj1 = axs[0,0].imshow(np.flipud(DEM_hs),cmap='Greys',alpha=1, interpolation='bilinear',extent=[DEM_xC[0,0],DEM_xC[0,-1],DEM_yC[0,0],DEM_yC[-1,0]])
demobj = axs[0,0].contourf(DEM_xC,DEM_yC,DEM_zC,cmap='gist_earth',alpha=0.5, levels=80)
axs[0,0].invert_yaxis()
axs[0,0].grid(c='k')
axs[0,0].set_title('a)')
for i in np.arange(n_stns):
    circ=plt.Circle((stn_xa[i],stn_ya[i]),10,color=colours(i),lw=0.5,ec='k')
    axs[0,0].add_artist(circ)

mmbeta=[np.floor(np.min(abs(incdpthout[:,:,6:8]))),np.ceil(np.max(abs(incdpthout[:,:,6:8])))]

# make beta vs depth plot
for i in np.arange(n_stns):
    axs[0,1].plot(GW_da,incdpthout[i,:,7],c=colours(i),linewidth=3)
axs[0,1].set_xlabel('Depth [m]')
axs[0,1].set_ylabel(r'$\beta$ [$\mu$Gal/m$_{H2O}]$')
axs[0,1].set_xlim([min(GW_da),max(GW_da)])
axs[0,1].set_ylim(mmbeta)
axs[0,1].grid()
axs[0,1].set_title('b)')

# make beta_z vs depth plot
for i in np.arange(n_stns):
    axs[1,0].plot(GW_da,-incdpthout[i,:,6],c=colours(i),linewidth=3)
axs[1,0].set_xlabel('Depth [m]')
axs[1,0].set_ylabel(r'$\beta_z$ [$\mu$Gal/m$_{H2O}]$')
axs[1,0].set_xlim([min(GW_da),max(GW_da)])
axs[1,0].set_ylim(mmbeta)
axs[1,0].grid()
axs[1,0].set_title('c)')

# make theta vs depth plot
theta=np.arccos(np.divide(-incdpthout[:,:,6],incdpthout[:,:,7]))*180/np.pi
for i in np.arange(n_stns):
    axs[1,1].plot(GW_da,theta[i,:],c=colours(i),linewidth=3)
axs[1,1].set_xlabel('Depth [m]')
axs[1,1].set_ylabel(r'$\theta$ [$^\circ$]')
axs[1,1].set_xlim([min(GW_da),max(GW_da)])
axs[1,1].grid()
axs[1,1].set_title('d)')
    
plt.savefig('Output/Rechy_betavsdepth1.pdf')

#%% profile plot
topoin=pd.read_csv('Other/topo_profile_fig_betavsdepth.csv')
plt.figure()
plt.plot(topoin['d'],topoin['z'],'k',lw=3)
plt.autoscale(enable=True, axis='both', tight=True)
plt.grid()
plt.xlabel('x-y distance [m]')
plt.ylabel('elevation [m]')

#%% effect of integration radius limit
incaresidout=[]
arl = np.logspace(-2,-.5,10)
dpth=2.5
for i in np.arange(n_stns):
#i=11
#if True:
    print('stn = '+str(i))
    incaresidout1=[]
    stn_x=stn_xa[i]
    stn_y=stn_ya[i]
    for ar in arl:
        out = Gravi4GW.Gravi4GW(DEM_path, stn_x, stn_y, dpth, accept_resid=ar, n_r=120, do_figs=False)
        incaresidout1.append(np.append(out,ar))
    incaresidout1 = np.array(incaresidout1).squeeze()
    incaresidout.append(incaresidout1)
incaresidout = np.array(incaresidout).squeeze()
#%% and plot it...
colours = cm.get_cmap('cividis_r',n_stns)
fig, ax = plt.subplots(nrows=1,ncols=1,sharex=False,sharey=False,figsize=(10,10))
mmbeta2=[np.floor(np.min(abs(incaresidout[:,:,6]))),np.ceil(np.max(abs(incaresidout[:,:,6])))]

# make beta_z vs depth plot
for i in np.arange(n_stns):
    ax.plot(incaresidout[i,:,8],-incaresidout[i,:,6],c=colours(i),linewidth=3)
ax.set_xlabel('$\epsilon$ [-]')
ax.set_ylabel(r'$\beta_z$ [$\mu$Gal/m$_{H2O}]$')
#ax.set_xlim([min(GW_da),max(GW_da)])
ax.set_ylim(mmbeta2)
ax.set_xscale('log')
ax.grid()
ax.invert_xaxis()
#ax.set_title('c)')

    
#plt.savefig('Output/Rechy_betavsacceptableresid1.pdf')


