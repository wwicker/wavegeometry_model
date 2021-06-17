# Python routines to prepare the basic state of a zonally averaged atmosphere
# based on the wavegeometry model by Nili Harnik and collegues
# May 2021, Wolfgang Wicker <wolfgang.wicker@env.ethz.ch>

import numpy as np
from scipy.ndimage import convolve1d
from scipy.interpolate import interp1d

## Define a grid class, methods and functions that use numpy.arrays for coordinates and dask.arrays for data
# EDIT: there is some problem with dask.array.apply_along_axis -> rewritten to use numpy.arrays for everything
# assume dimensions (time,z,y)

class Grid:
    def __init__(self,namelist,**kwargs):
        # define default parameters
        self.param = dict(ny=47,nz=71,nt=1,
                          top=15,bottom=0,wide=90,hem=1,
                          Ho=7000,om=7.27*10**(-5),ae=6.37*10**6,
                          Ro=287,g=9.81,kappa=0.2854,
                          alheight=10.5,alscale=2.5,alamp=1.5,
                          alwidth=15,alscaley=10,alampy=1,amrad=0,
                         )
        
        # read pairs of key and value from namelist
        # ignore lines that start with '#', use whitespaces as delimited
        f = open(namelist,'r')
        for line in f:
            if line.startswith('#'):
                continue
            else:
                l = line.strip().split(' ',1)
                try:
                    self.param.update({l[0]:int(l[1])})
                except ValueError:
                    self.param.update({l[0]:float(l[1])})
        f.close()
        
        # update parameters from kwargs
        self.param.update(kwargs)
        
        self.delz = (self.param['top'] - self.param['bottom'])/self.param['nz']*self.param['Ho']
        width = self.param['wide'] * np.pi / 180
        self.dely = width / (self.param['ny']-1)
        
        # y-coordinate (increasing latitude)
        if self.param['hem'] > 0: # Northern Hemisphere
            self.y = np.linspace(0,width,self.param['ny'])
        else:
            self.y = np.linspace(-np.pi/2,width-np.pi/2,self.param['ny'])
        # z-coordinate [bottom,top)
        self.z = np.arange(self.param['bottom']*self.param['Ho'],self.param['top']*self.param['Ho'],self.delz)
        # time-step coordinate
        self.t = np.arange(self.param['nt'])
        
        self.fcor = 2*self.param['om']*np.sin(self.y[np.newaxis,:])
        self.beta = 2*self.param['om']*np.cos(self.y[np.newaxis,:]) / self.param['ae']
        self.sn = np.sin(self.y[np.newaxis,:])
        self.cs = np.cos(self.y[np.newaxis,:])
        self.ez = np.exp(self.z[:,np.newaxis]/self.param['Ho'])
        
                
    def _fy(self,data):
        '''
            First y-derivative
        '''
        derv = lambda a: (a[:,:,2:] - a[:,:,:-2]) / 2 / self.dely
        result = derv(np.pad(data,((0,0),(0,0),(1,1)),mode='edge'))
        return result

    def _fyy(self,data):
        '''
            Second y-derivative
        '''
        derv = lambda a: (a[:,:,2:] + a[:,:,:-2] - 2*a[:,:,1:-1]) / self.dely / self.dely
        result = derv(np.pad(data,((0,0),(0,0),(1,1)),mode='edge'))
        return result
            
            
    def _fz(self,data):
        '''
            First z-derivative 
        '''
        derv = lambda a: (a[:,2:,:] - a[:,:-2,:]) / 2 / self.delz
        result = derv(np.pad(data,((0,0),(1,1),(0,0)),mode='edge'))
        return result
        
            
    def _fznz(self,data,N2):
        '''
            Second z-derivative taking into account density and N2 effects
            1/rho * (rho/N2*f_z)_z -> expand product rule in diff
        '''
        factor = self.ez / 2 / self.delz / self.delz
        roN2 = 1 / self.ez / N2
        
        diff = lambda a,roN2: ((roN2[:,2:,:]+roN2[:,1:-1,:])*a[:,2:,:] - 
                               (roN2[:,2:,:]+2*roN2[:,1:-1,:]+roN2[:,:-2,:])*a[:,1:-1,:] + 
                               (roN2[:,1:-1,:]+roN2[:,:-2,:])*a[:,:-2,:])
        result = diff(np.pad(data,((0,0),(1,1),(0,0)),mode='edge'),np.pad(roN2,((0,0),(1,1),(0,0)),mode='edge'))
        result = result * factor
        return result

            
            
            
def _interpolation(a,xi,x):
    '''
        Cubic spline interpolation
    '''
    f = interp1d(xi,a,kind='cubic',bounds_error=False,fill_value=np.nan)
    result = f(x)
    return result
    
def _2dinterpolation(a,yi,zi,y,z):
    '''
        That would require to load the 2d array into memory
    '''
    print('Not implemented yet.')
    return

def _filter(data,axis,kernel=np.array([0.25,0.5,0.25])):
    '''
        Smooth data along axis by convolution with kernel
    '''
    result = np.apply_along_axis(convolve1d,axis,data,kernel,mode='mirror')
                    
    return result
    
    
    
## User-level functions

def interp2z(grid,data,pressure,nfilter=None):
    '''
        Interpolate data from pressure-coordinates to vertical grid with constant log-pressure spacing
            
        pressure in Pa
    '''
    ps = np.max(pressure)
    zi = (-1)*np.log(pressure/ps)*grid.param['Ho']
    result = np.apply_along_axis(_interpolation,-2,data,zi,grid.z)
    
    if not(nfilter is None):
        for n in range(nfilter):
            result = _filter(result,-2)
    return result

    
def buoancy_frequency(grid,T):
    '''
        Calculate buoancy frequancy N2 from zonal mean temperature
    '''
    N2 = (grid._fz(T)+T*grid.param['kappa']/grid.param['Ho'])*grid.param['Ro']/grid.param['Ho']
    #for n in range(10):
    #    N2 = _filter(N2,-2)
    return N2
    
    
    
def pv_gradient(grid,U,N2):
    '''
        Calculate the PV gradient from zonal-mean zonal wind
        taking into account N2 effects for the vertical derivative
    '''
    #for n in range(2):
    #    U = _filter(U,-1)
    
    qy = U/grid.cs/grid.cs+grid._fy(U)*grid.sn/grid.cs
    # division by cosy produces out-of-bounds values
    qy[:,:,0] = 0
    qy[:,:,-1]  = 0
    qy = 1/grid.param['ae']**2 * (qy-grid._fyy(U)) - grid._fznz(U,N2)*grid.fcor*grid.fcor + grid.beta
    return qy
    
    
    
## PROCEDURE  

# Parse arguments: T-file, U-file, namelist text-file

# Read namelist and define model grid        

# Derive zonal-mean zonal wind and temperature from netcdf files
# - ensure that dimensions are (time,z,y)
# - data and model are expected to share the meridional grid
#   (regular grid spacing, one hemisphere, ny grid points including equator and pole, increasing latitudes)
#    -> this is for wide=90

# Change vertical coordinate to log pressure
# - Interpolate data to model grid
# - Apply 121 filter to smooth data in vertical (10 times or T, twice for U)

# Calculate buoancy frequency N2 from in-situ temperature

# Calculate PV gradient qy from zonal velocity
