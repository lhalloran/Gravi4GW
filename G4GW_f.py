# G4GW_f.py
# functions 
import numpy as np
import matplotlib.pyplot as plt

# define constants (all in SI units)
G=6.67408E-11 # gravitational constant
rho_H2O = 1000 # density of water

def quickplotterxyz(x,y,z):
    # quickplotterxyz(x,y,z)
    # Simple plotting of x,y,z arrays for debugging.
    # x,y,z : 2-D arrays of x, y, and elevation
    fig, axs = plt.subplots(nrows=1,ncols=3,sharex=True,sharey=True,figsize=(12,6)) 
    axs[0].imshow(x); axs[0].set_title('x, min='+str(int(min(x.flatten())))+', max='+str(int(max(x.flatten()))))
    axs[1].imshow(y); axs[1].set_title('y, min='+str(int(min(y.flatten())))+', max='+str(int(max(y.flatten()))))
    axs[2].imshow(z); axs[2].set_title('z, min='+str(int(min(z.flatten())))+', max='+str(int(max(z.flatten()))))

def point_maker(stn_x,stn_y,max_r,nr,dens_rad=8):
    # point_maker(stn_x,stn_y,max_r,nr,dens_rad=8)
    #
    # Returns x,y coordinates the element area represented by each point. Point distribution is 
    # radial with increasing distance from centre for each subsequent radius. Even radial spacing 
    # for each r.
    # 
    # Parameters:
    # stn_x, stn_y : coordinates of centre point
    # max_r : distance of furthest points from centre point
    # nr : number of radii (including r=0) for point creation
    # dens_rad : (optional) point density for first non-zero radius, must be multiple of 4
    #
    if dens_rad%4 != 0:
        print('#Gravi4GW: point_maker requires a radial density parameter dens_rad that is a multiple of 4')
        return 0
    rs=np.append(np.array(0.0),np.logspace(np.log10(0.1), np.log10(max_r), num=nr-1))
    log_fac=np.log10(rs[-1]/rs[-2])
    rsextra=rs[-1]*10**log_fac # for calculation of elemental areas of farthest points
    n=rs-rs
    n=dens_rad + 4*np.floor_divide(rs,10) #number of points for given r, increasing with radius
    n[0]=1 # one point at centre
    A = rs-rs # area represent by each point at given r
    for i in np.arange(1,nr-1):
        r_out = (rs[i+1]+rs[i])/2 # outer radius of element
        r_in = (rs[i]+rs[i-1])/2 # inter radius of element
        A[i] = np.pi*(r_out**2-r_in**2)/n[i]
    A[0] = np.pi*((rs[0]+rs[1])/2)**2
    A[-1]=np.pi*((rsextra+rs[-1])**2-(rs[-1]+rs[-2])**2)/(4*n[-1])
    # make r,phi coordinates, referenced to centre point
    pt_r=np.array([])
    pt_phi=np.array([])
    pt_A=np.array([])
    for i in np.arange(nr):
        for j in np.arange(n[i]):
            pt_r = np.append(pt_r,rs[i])
            #print(i)
            #print(pt_r)
            pt_phi = np.append(pt_phi,2*np.pi*j/n[i])
            pt_A=np.append(pt_A,A[i])
    # convert to x,y cords:
    x=stn_x+np.multiply(pt_r,np.cos(pt_phi))
    y=stn_y+np.multiply(pt_r,np.sin(pt_phi))
    return x,y,pt_A

def dg(xyz_stn,xyz_pt,dm,debug=False):
    # dg(xyz_stn,xyz_pt,dm,debug=False) is a function that returns the delta gravity vector for a station point and 
    # a grid point, given a change in mass of dm at the point
    #
    # xyz_stn is local coordiates (in m) of gravity station. numpy array of size 3.
    # xyz_pt is local coordiates (in m) of point on GW table grid. numpy array of size 3.
    # dm is (near infinitessimal) change in mass of water at the point on the grid
    # (optional) debug=True will print internal values for each call of the function
    dxyz = xyz_pt-xyz_stn;
    dx=dxyz[0]
    dy=dxyz[1]
    dz=dxyz[2]
    d = np.sqrt(dx**2+dy**2+dz**2) #distance between station and point
    dg = G*dm/d**2
    theta = np.arctan2(dz,np.sqrt(dx**2+dy**2)) # angle below xy surface
    phi = np.arctan2(dy,dx)    
    dgx=dg*np.cos(theta)*np.sin(phi)
    dgy=dg*np.cos(theta)*np.cos(phi)
    dgz=dg*np.sin(theta)
    dgxyz=np.array([dgx,dgy,dgz])
    if(debug):
        print(dxyz)
        print([d,dg,theta,phi])
        print(dgxyz)
    return dgxyz
############################################################################