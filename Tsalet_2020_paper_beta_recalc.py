import Gravi4GW
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm
from osgeo import gdal

DEM_path = 'DEMs/RechyDEM_swisstopo_2m.tif'
stns_in = pd.read_csv('Other/Tsalet_HP_Paper_Gravi_Stns.csv')
nstns = stns_in.shape[0]
gravimeter_sensor_height = 0.48
GW_ds = np.arange(2.5,10.1,1.25) + gravimeter_sensor_height # +0.48 is height of sensor above ground
ndepths=np.size(GW_ds)

file_out = 'Output/Tsalet_2020_paper_beta_recalc_data_out.csv'
#%%
data_out = []
for sn,stn_x,stn_y,dgog in zip(stns_in['Station'],stns_in['X'],stns_in['Y'],stns_in['dg 2019']):
    for GW_d in GW_ds: 
        print(str(stn_x)+', '+str(stn_y))
        output = Gravi4GW.Gravi4GW(DEM_path, stn_x, stn_y, GW_d, accept_resid=0.02, n_r=80, do_figs=False)
        output = np.append([sn,dgog],output)
        data_out.append(output)
data_out = np.array(data_out).squeeze()
#%%
nrows=np.shape(data_out)[0]
extracalcs=np.zeros((nrows,3))
extracalcs[:,0] = np.arccos(np.divide(-data_out[:,8],data_out[:,9]))*180/np.pi
extracalcs[:,1]=-np.divide(data_out[:,1],data_out[:,8])
extracalcs[:,2]=data_out[:,1]/41.93
data_out2 = np.append(data_out,extracalcs,axis=1)

dfcols = ['Station', '2019 delta g', 'x [m]','y [m]','z [m]','eff. depth [m]','beta (x-component) [uGal/m]','beta (y-component) [uGal/m]','beta (z-component) [uGal/m]','beta [uGal/m]','theta_beta','delta H (using new beta)', 'delta H (beta=41.93)']
outputdf = pd.DataFrame(data=data_out2,columns=dfcols)
outputdf.to_csv(file_out,index=False)
#%%
plt.figure(figsize=(6,8))
colours = cm.get_cmap('nipy_spectral',nstns)
for i in np.arange(0,nrows,ndepths):
    cnow=np.array(colours(int(i/ndepths)))*0.86; cnow[3]=1
    plt.plot(outputdf['eff. depth [m]'][i:i+ndepths]-gravimeter_sensor_height, outputdf['delta H (using new beta)'][i:i+ndepths], c=cnow, lw=2)
    plt.plot(0,outputdf['delta H (beta=41.93)'][i],'>',c=cnow)
    if outputdf['Station'][i]==9:
        textxpos=0.3
    elif outputdf['Station'][i]==11:
        textxpos=1.4
    elif outputdf['Station'][i]==8:
        textxpos=0.5
    elif outputdf['Station'][i]==10:
        textxpos=1.1
    else:
        textxpos=0.8
    plt.text(textxpos,outputdf['delta H (beta=41.93)'][i],'G'+str(int(outputdf['Station'][i])),c=cnow,horizontalalignment='center',verticalalignment='center')
plt.xlabel(r'assumed h$_{eff}$ [m]',fontsize=14)
plt.ylabel(r'$\Delta h$ [m$_{H2O}$]',fontsize=14)
plt.xticks(ticks=[0,2.5,5,7.5,10],labels=['BPA','2.5','5','7.5','10'])
plt.grid()

#%% map with points
DEM_in = gdal.Open(DEM_path, gdal.GA_ReadOnly) 
DEM_z = np.array(np.float64(DEM_in.ReadAsArray()))
DEM_geotransform = DEM_in.GetGeoTransform()
xy_inds = np.indices((DEM_in.RasterXSize, DEM_in.RasterYSize))
DEM_x = DEM_geotransform[0] + DEM_geotransform[1]*xy_inds[0] + DEM_geotransform[2]*xy_inds[1]
DEM_y = DEM_geotransform[3] + DEM_geotransform[4]*xy_inds[0] + DEM_geotransform[5]*xy_inds[1]
DEM_x=DEM_x.transpose()
DEM_y=DEM_y.transpose()
mmx=np.array([2605957.0+120, 2606937.0])
mmy=np.array([1115619.0+120, 1116599.0])+20
xinds = np.where(np.logical_and(np.min(DEM_x,axis=0) > mmx[0], np.max(DEM_x,axis=0) < mmx[1]))[0]
yinds = np.where(np.logical_and(np.min(DEM_y,axis=1) > mmy[0], np.max(DEM_y,axis=1) < mmy[1]))[0]
DEM_xC = DEM_x[yinds,:][:,xinds]
DEM_yC = DEM_y[yinds,:][:,xinds]
DEM_zC = DEM_z[yinds,:][:,xinds]

fig, axs = plt.subplots(nrows=1,ncols=1,sharex=False,sharey=False,figsize=(6,6))
DEM_hs = Gravi4GW.hillshade(DEM_zC,45,20)
cbobj1 = axs.imshow(np.flipud(DEM_hs),cmap='Greys',alpha=1, interpolation='bilinear',extent=[DEM_xC[0,0],DEM_xC[0,-1],DEM_yC[0,0],DEM_yC[-1,0]])
demobj = axs.contourf(DEM_xC,DEM_yC,DEM_zC,cmap='gist_earth',alpha=0.5, levels=80)
axs.invert_yaxis()
axs.grid(c='k')

#for sn,stn_x,stn_y,dgog in zip(stns_in['Station'],stns_in['X'],stns_in['Y'],stns_in['dg 2019']):
dxtext=0
dytext=-15
for i in np.arange(0,nstns):
    cnow=np.array(colours(i))*0.86; cnow[3]=1
    circ=plt.Circle((stns_in['X'][i],stns_in['Y'][i]),10,color=colours(i),lw=0.5,ec='k')
    axs.add_artist(circ)
    betaznow = np.array(outputdf[outputdf['eff. depth [m]']==5.48]['beta (z-component) [uGal/m]'])[i]
    strnow='G'+str(stns_in['Station'][i]) + '|' + "{:.1f}".format(-betaznow)
    axs.text(stns_in['X'][i]+dxtext,stns_in['Y'][i]+dytext,strnow,c='w',horizontalalignment='center',verticalalignment='center')
axs.text(mmx[0]+75,mmy[0]+20,r'Stn|$\beta_z$',c='w',horizontalalignment='center',verticalalignment='center')