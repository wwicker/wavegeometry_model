# Python routines to solve for quasi-geostrophic wave solutions in a given basic state
# spherical coordinate using Plumb's formulation (streamfunction=geopotential height/coriolis force)
# based on the wavegeometry model by Nili Harnik and collegues
# May 2021, Wolfgang Wicker <wolfgang.wicker@env.ethz.ch>

import numpy as np
import scipy.sparse as sparse
import scipy.sparse.linalg as linalg



# The model equation is given equation (2) in Harnik & Lindzen (2001) but adapted for sphere
# instead of a beta plane similiar equation (C1) ibid.
#   A * phi = f
# At the bottom phi is given by the forcing f, at the top and at the equatorial boundary, wave reflection is prohibited by 
# linear Rayleigh friction and Newtonian cooling.

# Using second order-accurate finite differences to discretize that equation derives in a sparse matrix equation
# with five non-zero diagononals. For solving that equation for a 3D basic state (time,z,y) it is looped
# over the time-dimension. To obtain 1D vector phi and a 2D matrix A the z- and y-dimension are stacked
# with y changing fastest. The 2D sparse matrix can be solved using scipy.sparse.linalg.

# Both the functions to define the sponge layer and to compute the diagnostics are based on the Grid class in basicstate.py.
# The involved methods to compute derivatives use second order-accurate finite differences. That leave room for improvement.



### USER-LEVEL FUNCTIONS to define the sponge layer        
        
def rayd(grid):
    # set z-dependence of linear damping coefficient (damping at the model top)
    alheight = grid.param['alheight'] * grid.param['Ho']
    alscale = grid.param['alscale'] * grid.param['Ho']
    zrayd = np.tanh((grid.z-alheight) / alscale) - np.tanh((grid.z[0]-alheight) / alscale)
    zrayd = 0.5 * grid.param['alamp']/86400 * zrayd
    
    #set y-dependence of linear damping coefficient (damping at the equator)
    if grid.param['hem'] > 0: # Northern Hemisphere
        alwidth = np.abs(grid.param['alwidth']) * np.pi / 180
        alscaley = np.abs(grid.param['alscaley']) * np.pi / 180
        yrayd = np.tanh((grid.y[-1]-alwidth) / alscaley) - np.tanh((grid.y-alwidth) / alscaley)
    else:
        alwidth = - np.abs(grid.param['alwidth']) * np.pi / 180
        alscaley = - np.abs(grid.param['alscaley']) * np.pi / 180
        yrayd = np.tanh((grid.y-alwidth) / alscaley) - np.tanh((grid.y[0]-alwidth) / alscaley)
    yrayd = 0.5 * grid.param['alampy']/86400 * yrayd
    
    coeff = zrayd[:,np.newaxis] + yrayd
    return coeff[:,:]
    
    
def alpha(grid):
    r = rayd(grid)
    
    # I don't understand the z-dependent factor
    coeff = r + grid.param['amrad'] * (0.45*np.exp(-(grid.z-50000)**2/13000**2)+0.05)[:,np.newaxis]
    return coeff
        
        
def damp(grid,U,s=2,gammac=0.05,freq=0):
    # get complex part of phase speed
    gamma = 1j * gammac * grid.param['ae'] * grid.cs / s /86400
    # get real part of phase speed from frequency which is given in days**(-1)
    cph = freq * grid.param['ae'] * grid.cs / s /86400
        
    coeff =  1j * s * (U-cph-gamma) / grid.param['ae'] / grid.cs
    return coeff
            
              
        
        
class Matrix:
    def __init__(self,grid,s=2):
        
        self.s = s
        
        # Define grid attributes
        self.ae = grid.param['ae']
        self.nt = grid.param['nt']
        self.nz = grid.param['nz']
        self.ny = grid.param['ny']
        
        self.delz = grid.delz
        self.dely = grid.dely
        
        self.fcor = grid.fcor
        self.sn = grid.sn
        self.cs = grid.cs
        self.ez = grid.ez
        
    
    ### INTERIOR COEFFICIENTS
    # need to set boundary coefficients to zero
    
    def _a5(self,qy,roN2,damp,rayd,alpha):
        '''
            Main diagonal
        '''
        coeff = 1j*self.s/self.ae/self.cs[:,1:-1]*qy[:,1:-1,1:-1]
        coeff = coeff - self.s**2/self.ae**2/self.cs[:,1:-1]**2 * (damp[:,1:-1,1:-1]+rayd[1:-1,1:-1])
        coeff = coeff - ((roN2[:,2:,1:-1]+2*roN2[:,1:-1,1:-1]+roN2[:,:-2,1:-1])/2/self.delz**2 * 
                         self.fcor[:,1:-1]**2*self.ez[1:-1,:]*(damp[:,1:-1,1:-1]+alpha[1:-1,1:-1]))
        coeff = coeff - 2*(damp[:,1:-1,1:-1]+rayd[1:-1,1:-1])/self.dely**2/self.ae**2
        return np.pad(coeff,((0,0),(1,1),(1,1)),'constant',constant_values=(1,))
       
        
    def _a4(self,damp,rayd):
        '''
            First diagonal below the main diagonal
        '''
        coeff = (damp[:,1:-1,1:-1]+rayd[1:-1,1:-1])/self.dely**2/self.ae**2
        coeff = coeff + (damp[:,1:-1,1:-1]+rayd[1:-1,1:-1])*self.sn[:,1:-1]/self.cs[:,1:-1]/self.ae**2/2/self.dely
        coeff = coeff - (rayd[1:-1,2:]-rayd[1:-1,:-2])/4/self.dely**2/self.ae**2
        return np.pad(coeff,((0,0),(1,1),(1,1)),'constant',constant_values=(0,))
    
    
    def _a6(self,damp,rayd):
        '''
            First diagonal above the main diagonal
        '''
        coeff = (damp[:,1:-1,1:-1]+rayd[1:-1,1:-1])/self.dely**2/self.ae**2
        coeff = coeff - (damp[:,1:-1,1:-1]+rayd[1:-1,1:-1])*self.sn[:,1:-1]/self.cs[:,1:-1]/self.ae**2/2/self.dely
        coeff = coeff + (rayd[1:-1,2:]-rayd[1:-1,:-2])/4/self.dely**2/self.ae**2
        return np.pad(coeff,((0,0),(1,1),(1,1)),'constant',constant_values=(0,))
    
    
    def _a2(self,N2,roN2,damp,alpha):
        '''
            nyth diagonal below the main diagonal
        '''
        coeff = (roN2[:,1:-1,1:-1]+roN2[:,:-2,1:-1])/2/self.delz**2 * (self.fcor[:,1:-1]**2 * self.ez[1:-1,:] *
                                                                        (damp[:,1:-1,1:-1]+alpha[1:-1,1:-1]))
        coeff = coeff - self.fcor[:,1:-1]**2/N2[:,1:-1,1:-1] * (alpha[2:,1:-1]-alpha[:-2,1:-1])/4/self.delz**2
        return np.pad(coeff,((0,0),(1,1),(1,1)),'constant',constant_values=(0,))
        
        
    def _a8(self,N2,roN2,damp,alpha):
        '''
            nyth diagonal above the main diagonal
        '''
        coeff = (roN2[:,1:-1,1:-1]+roN2[:,2:,1:-1])/2/self.delz**2 * (self.fcor[:,1:-1]**2 * self.ez[1:-1,:] *
                                                                        (damp[:,1:-1,1:-1]+alpha[1:-1,1:-1]))
        coeff = coeff + self.fcor[:,1:-1]**2/N2[:,1:-1,1:-1] * (alpha[2:,1:-1]-alpha[:-2,1:-1])/4/self.delz**2
        return np.pad(coeff,((0,0),(1,1),(1,1)),'constant',constant_values=(0,))
        
        
        
    def _stack(self,data):
        '''
            stack z- and y-dimension to (z_0,y_0),(z_0,y_1), ..., (z_nz,y_ny-1),(z_nz,y_ny)
            
            data should have dimensions (nt,nz,ny)
        '''
        return np.reshape(data,(self.nt,-1),order='C')
    
    
    def _unstack(self,data):
        '''
            reshape to dimensions (nt,nz,ny)
        '''
        return np.reshape(data,(self.nt,self.nz,self.ny),order='C')
    
    
    ### USER-LEVEL METHODS to define forcing, compute matrix coefficients, and solve the matrix equation
    
    def A(self,qy,N2,damp,rayd,alpha):
        '''   
            matrix A in the QG linear system
            
            CHECK THE MATRIX EQUATION
            Not sure whether I have to cut the first or the last elements
        '''
        roN2 = 1 / self.ez / N2
        diagonals = [self._stack(self._a2(N2,roN2,damp,alpha))[:,self.ny:],
                     self._stack(self._a4(damp,rayd))[:,1:],
                     self._stack(self._a5(qy,roN2,damp,rayd,alpha))[:,:],
                     self._stack(self._a6(damp,rayd))[:,:-1],
                     self._stack(self._a8(N2,roN2,damp,alpha))[:,:-self.ny],]
        offsets = [-self.ny,-1,0,1,self.ny]
        return diagonals, offsets
    
    
    def solve(self,A,f):
        '''
            define sparse matrix from list of diagonals A[0] and offsets A[1]
            solve matrix equation A*phi=f
        '''
        phi = []
        for i in range(self.nt):
            data = [a[i,:] for a in A[0]]
            diag = sparse.diags(data,A[1],format='csr')
            phi.append(linalg.spsolve(diag,f))
        phi = np.stack(phi,axis=0)
        phi = self._unstack(phi)
        return phi
    
    
    def forcing(self):
        '''
            with the boundary values of _a5 set to 1, the forcing defines phi at the level z_0
        '''
        # Forcing is zero everywhere apart from the lowest level
        f = np.zeros((self.nz,self.ny))
        # Meridional structure of geopotential height is set to constant 1
        # divide by coriolis to obtain streamfunction
        f[0,1:-1] = 1 / self.fcor[:,1:-1]
        return f.flatten()
    

    
## User level functions to compute Diagnostics    
    
    
def FN2(grid,N2):
    '''
        contribution of density and N2 to the index of refraction
    '''
    f = np.sqrt(grid.ez*N2)
    result = 1 / N2 / np.sqrt(grid.ez) * grid._fznz(f,N2)
    return result
    
    
def l2(grid,phi):
    '''
        Meridional wavenumber squared
    '''
    result = - np.real((grid._fyy(phi)-grid.sn/grid.cs*grid._fy(phi))/phi)
    result[np.isinf(result)] = 0
    result[np.isnan(result)] = 0
    return result
        
        
def m2(grid,phi,N2):
    '''
        Vertical wavenumber squared
    '''
    psizz = np.sqrt(N2/grid.ez) * (grid._fznz(phi,N2)-FN2(grid,N2)*phi)
    result = - np.real(psizz/phi) * np.sqrt(grid.ez * N2)
    result[np.isinf(result)] = 0
    result[np.isnan(result)] = 0
    return result
    
    
def n2ref(grid,phi,N2):
    '''
        Refractive index squared
    '''
    result = m2(phi,N2) + N2 / grid.fcor**2 * l2(phi)
    result[np.isnan(result)] = 0
    return result
    
    
def n2ref_alt(grid,N2,qy,U,s,freq):
    '''
        Alternative formula to calculate the refractive index squared
    '''
    cph = freq * grid.param['ae'] * grid.cs / s/86400
    result = qy / (U-cph) * grid.param['ae']**2
    result = result - s**2 / grid.cs**2
    result = result + FN2(grid,N2) * grid.param['ae']**2 * grid.fcor**2
    result = N2 * result
    result[np.isinf(result)] = 0
    result[np.isnan(result)] = 0
    return result
        
        
def amplitude(grid,phi):
    '''
        DO I HAVE TO DIVIDE BY FCOR to get amplitude in units of geopotential height?
    '''
    result = np.abs(phi)
    return result
    
    
def phase(grid,phi):
    '''
        Phase angle 
    '''
    result = np.arctan2(np.imag(phi),np.real(phi))
    return result
        
        
def EPflux(grid,phi,N2,s):
    '''
        Meridional and vertical component of the quasigeostrophic Eliassen-Palm Flux
    '''
    epfy = np.imag(np.conjugate(phi)*grid._fy(phi)) / grid.ez / 2 / grid.param['ae'] * s
    epfz = np.imag(np.conjugate(phi)*grid._fz(phi)) / grid.ez / N2 * grid.fcor**2  / 2 * s
    return epfy, epfz
    
    
def enstrophy(grid,phi,N2,s):
    '''
        Squared potential vorticity perturbation
    '''
    phixx = - phi * s**2 / grid.param['ae']**2 / grid.cs**2
    qprime = grid.fcor**2 * grid._fznz(phi,N2) + phixx 
    qprime = qprime + grid._fyy(phi) / grid.param['ae']**2 + grid._fy(phi) * grid.sn / grid.cs / grid.param['ae']**2
    result = np.abs(qprime)**2
    return result
    
    


## PROCEDURE 

# initialize model instance

# read basic state

# define sponge field similiar to the basic state fields

# stack fields

# define and solve the 2D matrix equation

# compute diagnostics
#  - zonal and meridional wavenumber using the second derivative of phi
#  - EP flux depending on N2, phi, and the first derivative of phi
#  - refractive index depending of qy, U, cph, N2
#    (refractive index can also be computed from m2, l2 -> from phi)

#  - wave amplitude and phase from phi

#  - enstrophy & wave activity density depending on phi
#    (remember to transform streamfunction back to geopotential height)