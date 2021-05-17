# Python routines to prepare the basic state of a zonally averaged atmosphere
# based on the wavegeometry model by Nili Harnik
# May 2021, Wolfgang Wicker <wolfgang.wicker@env.ethz.ch>

import numpy as np

## Define a grid class, methods and functions that use numpy.arrays for coordinates and dask.arrays for data

class grid:
    def __init__(self,namelist,**kwargs):
        # define default parameters
        self.param = dict(ny=47,nz=71,nt=1
                          top=15,bottom=0,wide=90,hem=1,
                          Ho=7000,om=7.27*10**(-5),ae=6.37*10**6)
        
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
        
        
        def interp2z(self,data,pressure):
            '''
                Interpolate data from pressure-coordinates to vertical grid with constant log-pressure spacing
            '''
            
        def partial_y(self,data):
            '''
                First y-derivative
            '''

        def partial_yy(self,data):
            '''
                Second second y-derivative
            '''
            
            
        def partial_z(self,data):
            '''
                First z-derivative 
            '''
            
        def partial_zz(self,data,N2):
            '''
                Second z-derivative taking into account density and N2 effects
            '''

def spline_interpolation()            

def convolution_filter(data,axis,kernel=np.array[0.25,0.5,0.25]):
    '''
        Smooth data along axis by convolution with kerne
        
        - first and last element should not be affected (repeat those elements before convolution)
    '''
                                

def buoancy_frequency(T,grid):
    '''
        Calculate buoancy frequancy N2 from zonal mean temperature
    '''
    
    
    
def PV_grad(U,N2,grid):
    '''
        Calculate the PV gradient from zonal-mean zonal wind
        taking into account N2 effects for the vertical derivative
    '''
    
    
    
## PROCEDURE  

# Parse arguments: T-file, U-file, namelist text-file

# Read namelist and define model grid        

# Derive zonal-mean zonal wind and temperature from netcdf files

# Change vertical coordinate to log pressure

# Interpolate data to model grid

# Apply 121 filter to smooth data in vertical and meridional dimension

# Define routines for vertical and meridional derivatives

# Calculate buoancy frequency N2 from in-situ temperature

# Calculate PV gradient qy from zonal velocity