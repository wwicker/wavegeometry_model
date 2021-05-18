# Python routines to prepare the basic state of a zonally averaged atmosphere
# based on the wavegeometry model by Nili Harnik
# May 2021, Wolfgang Wicker <wolfgang.wicker@env.ethz.ch>

import numpy as np
import dask.array as da
from scipy.ndimage import convolve1d
from scipy.interpolate import interp1d

## Define a grid class, methods and functions that use numpy.arrays for coordinates and dask.arrays for data
# assume dimensions (time,z,y)

class grid:
    def __init__(self,namelist,**kwargs):
        # define default parameters
        self.param = dict(ny=47,nz=71,nt=1,
                          top=15,bottom=0,wide=90,hem=1,
                          Ho=7000,om=7.27*10**(-5),ae=6.37*10**6,
                          Ro=287,g=9.81)
        
        # read pairs of key and value from namelist
        # ignore lines that start with '#', use whitespaces as delimited
        f = open(namelist,'r')
        for line in f:
            if line.startswith('#'):
                continue
            else:
                l = line.strip().split(' ',1)
                self.param.update({l[0]:int(l[1])})
        f.close
        
        # update parameters from kwargs
        self.param.update(kwargs)
        
        self.delz = (self.param['top'] - self.param['bottom'])/self.param['nz']*self.param['Ho']
        width = self.param['wide'] * np.pi / 180
        self.dely = width / (self.param['ny']-1)
        
        # y-coordinate (increasing latitude)
        if hem > 0: # Northern Hemisphere
            self.y = np.linspace(0,width,self.param['ny'])
        else:
            self.y = np.linspace(-np.pi/2,width-np.pi/2,self.param['ny'])
        # z-coordinate [bottom,top)
        self.z = np.arange(self.param['bottom']*self.param['Ho'],self.param['top']*self.param['Ho'],delz)
        # time-step coordinate
        self.t = np.arange(self.param['nt'])
        
        self.fcor = 2*self.param['om']*np.sin(self.y)
        self.beta = 2*self.param['om']*np.cos(self.y) / self.param['ae']
        
        
        def interp2z(self,data,pressure,n_filter=None):
            '''
                Interpolate data from pressure-coordinates to vertical grid with constant log-pressure spacing
                
                pressure in Pa, assume 
            '''
            ps = np.max(pressure)
            zi = np.log(-pressure/ps)*self.param['Ho']
            result = data.apply_along_axis(interpolation,zi,self.z)
            if not(nfilter is None):
                for n in range(n_filter):
                    result = _filter(result,-2)
            return result
                
            
        def partial_y(self,data):
            '''
                First y-derivative
            '''
            derv = lambda a: (a[:,:,2:] - a[:,:,:-2]) / 2 / self.dely
            result = data.map_overlap(derv,depth=1,boundary='reflect')
            return result

        def partial_yy(self,data):
            '''
                Second second y-derivative
            '''
            derv = lambda a: (a[:,:,2:] + a[:,:,:-2] - 2*a[:,:,1:-1]) / self.dely / self.dely
            result = data.map_overlap(derv,depth=1,boundary='reflect')
            
            
        def partial_z(self,data):
            '''
                First z-derivative 
            '''
            derv = lambda a: (a[:,2:,:] - a[:,:-2,:]) / 2 / self.delz
            result = data.map_overlap(derv,depth=1,boundary='nearest')
            return result
        
            
        def partial_zz(self,data,N2):
            '''
                Second z-derivative taking into account density and N2 effects
                (rho/N2*f_z)_z -> expand product rule in diff
            '''
            factor = np.exp(self.z/self.param['Ho']) / 2 / self.delz / self.delz
            rho = np.exp(-self.z/self.param['Ho'])
            roN2 = 1 / N2 * rho[:,np.newaxis]
            
            diff = lambda a,roN2: ((roN2[:,2:,:]+roN2[:,1:-1,:])*a[:,2:,:] - 
                                   (roN2[:,2:,:]+2*roN2[:,1:-1,:]+roN2[:,:-2,:])*a[:,1:-1,:] + 
                                   (roN2[:,1:-1,:]+roN2[:,:-2,:])*a[:,:-2,:])
            result = data.map_overlap(diff,depth=1,boundary=nearest)
            result = result * factor[:,np.newaxis]
            return result

            
            
            
def _interpolation(a,xi,x):
    f = interp1d(xi,a,kind='cubic')
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
    result = data.apply_along_axis(convolve1d,axis,kernel,mode='mirror')
    return result
    
    
    
## User-level functions             

def buoancy_frequency(T,grid):
    '''
        Calculate buoancy frequancy N2 from zonal mean temperature
    '''
    N2 = (grid.partial_z(T)+T*2/7/grid.param['Ho'])*Ro/Ho
    for n in range(10):
        N2 = _filter(N2,-2)
    return N2
    
    
    
def PV_grad(U,N2,grid):
    '''
        Calculate the PV gradient from zonal-mean zonal wind
        taking into account N2 effects for the vertical derivative
    '''
    for n in range(2):
        U = _filter(U,-1)
    tany = np.tan(grid.y)
    cosy = np.cos(grid.y)
    
    qy = U/cosy/cosy+grid.partial_y(U)*tany+grid.partial_yy(U)-grid.partial_zz(U)*grid.fcor*grid.fcor
    qy = 1/grid.param['ae']**2 * qy + grid.beta
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