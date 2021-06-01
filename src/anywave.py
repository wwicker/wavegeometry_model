# Python routines to solve for quasi-geostrophic wave solutions in a given basic state
# spherical coordinate using Plumb's formulation (streamfunction=geopotential height/coriolis force)
# based on the wavegeometry model by Nili Harnik and collegues
# May 2021, Wolfgang Wicker <wolfgang.wicker@env.ethz.ch>

import numpy as np
import scipy.sparse as sparse
import scipy.sparse.linalg as linalg


# Define a Sponge class for all damping coefficients, define a Matrix class for the matrix equation
# all damping and matrix coefficients are computed as 3D array (time,z,y)
# stack z- and y-dimension, add all coefficients and solve the maxtrix equation A*phi=forcing
# use scipy.sparse.linalg to solve the equation for the sparce matrix A
# the Matrix for second order accurate centered differences contains 5 diagonals

class Sponge:
    def __init__(self,namelist,**kwargs):
        # define default parameters
        self.param = dict(ny=47,nz=71,nt=1,
                          top=15,bottom=0,wide=90,hem=1,
                          Ho=7000,ae=6.37*10**6,
                          alheight=10.5,alscale=2.5,alamp=1.5,
                          alwidth=15,alscaley=10,alampy=1,
                          amrad=0,freq=0,gammac=0.05,s=2)
        
        # read pairs of key and value from namelist
        # ignore lines that start with '#', use whitespaces as delimited
        f = open(namelist,'r')
        for line in f:
            if line.startswith('#'):
                continue
            else:
                l = line.strip().split(' ',1)
                self.param.update({l[0]:int(l[1])})
        f.close()
        
        # update parameters from kwargs
        self.param.update(kwargs)
        
        delz = (self.param['top'] - self.param['bottom'])/self.param['nz']*self.param['Ho']
        width = self.param['wide'] * np.pi / 180
        # y-coordinate (increasing latitude)
        if self.param['hem'] > 0: # Northern Hemisphere
            self.y = np.linspace(0,width,self.param['ny'])
        else:
            self.y = np.linspace(-np.pi/2,width-np.pi/2,self.param['ny'])
        # z-coordinate [bottom,top)
        self.z = np.arange(self.param['bottom']*self.param['Ho'],self.param['top']*self.param['Ho'],delz)
        
        
    def rayd(self):
        # set z-dependence of linear damping coefficient (damping at the model top)
        alheight = self.param['alheight'] * self.param['Ho']
        alscale = self.param['alscale'] * self.param['Ho']
        zrayd = np.tanh((self.z-alheight) / alscale) - np.tanh((self.z[0]-alheight) / alscale)
        zrayd = 0.5 * self.param['alamp']/86400 * zrayd
        
        #set y-dependence of linear damping coefficient (damping at the equator)
        if self.param['hem'] > 0: # Northern Hemisphere
            alwidth = np.abs(self.param['alwidth']) * np.pi / 180
            alscaley = np.abs(self.param['alscaley']) * np.pi / 180
            yrayd = np.tanh((self.y[-1]-alwidth) / alscaley) - np.tanh((self.y-alwidth) / alscaley)
        else:
            alwidth = - np.abs(self.param['alwidth']) * np.pi / 180
            alscaley = - np.abs(self.param['alscaley']) * np.pi / 180
            yrayd = np.tanh((self.y-alwidth) / alscaley) - np.tanh((self.y[0]-alwidth) / alscaley)
        yrayd = 0.5 * self.param['alampy']/86400 * yrayd
        
        coeff = zrayd[:,np.newaxis] + yrayd
        return coeff[:,:]
    
    
    def alpha(self):
        rayd = self.rayd()

        # I don't understand the z-dependent factor
        coeff = rayd + self.param['amrad'] * (0.45*np.exp(-(self.z-50000)**2/13000**2)+0.05)[:,np.newaxis]
        return coeff
        
        
    def damp(self,U):
        # get complex part of phase speed
        gamma = 1j * self.param['gammac'] * self.param['ae'] * np.cos(self.y) / self.param['s']/86400
        # get real part of phase speed from frequency which is given in days**(-1)
        cph = self.param['freq'] * self.param['ae'] * np.cos(self.y) / self.param['s']/86400
        
        coeff =  1j * self.param['s'] * (U-cph-gamma) / self.param['ae'] / np.cos(self.y)
        return coeff
            
              
        
        
class Matrix:
    def __init__(self,namelist,**kwargs):
        # define default parameters
        self.param = dict(ny=47,nz=71,nt=1,
                          top=15,bottom=0,wide=90,hem=1,
                          Ho=7000,om=7.27*10**(-5),ae=6.37*10**6,
                          s=2)
        
        # read pairs of key and value from namelist
        # ignore lines that start with '#', use whitespaces as delimited
        f = open(namelist,'r')
        for line in f:
            if line.startswith('#'):
                continue
            else:
                l = line.strip().split(' ',1)
                self.param.update({l[0]:int(l[1])})
        f.close()
        
        # update parameters from kwargs
        self.param.update(kwargs)
        
        self.delz = (self.param['top'] - self.param['bottom'])/self.param['nz']*self.param['Ho']
        width = self.param['wide'] * np.pi / 180
        self.dely = width / (self.param['ny']-1)
        
        # y-coordinate (increasing latitude)
        if self.param['hem'] > 0: # Northern Hemisphere
            y = np.linspace(0,width,self.param['ny'])
        else:
            y = np.linspace(-np.pi/2,width-np.pi/2,self.param['ny'])
        # z-coordinate [bottom,top)
        z = np.arange(self.param['bottom']*self.param['Ho'],self.param['top']*self.param['Ho'],self.delz)
        
        self.fcor = 2*self.param['om']*np.sin(y[np.newaxis,:])
        self.sn = np.sin(y[np.newaxis,:])
        self.cs = np.cos(y[np.newaxis,:])
        self.ez = np.exp(z[:,np.newaxis]/self.param['Ho'])
        

    
    
    ## INTERIOR COEFFICIENTS
    # need to set boundary coefficients to zero
    def _a5(self,qy,roN2,damp,rayd,alpha):
        '''
            Hauptdiagonale
        '''
        coeff = 1j*self.param['s']/self.param['ae']/self.cs[:,1:-1]*qy[:,1:-1,1:-1]
        coeff = coeff - self.param['s']**2/self.param['ae']**2/self.cs[:,1:-1] * (damp[:,1:-1,1:-1]+rayd[1:-1,1:-1])
        coeff = coeff - ((roN2[:,2:,1:-1]+2*roN2[:,1:-1,1:-1]+roN2[:,:-2,1:-1])/2/self.delz**2 * 
                         self.fcor[:,1:-1]**2*self.ez[1:-1,:]*(damp[:,1:-1,1:-1]+alpha[1:-1,1:-1]))
        coeff = coeff - 2*(damp[:,1:-1,1:-1]+rayd[1:-1,1:-1])/self.dely**2/self.param['ae']**2
        return np.pad(coeff,((0,0),(1,1),(1,1)),'constant',constant_values=(1,))
       
        
    def _a4(self,damp,rayd):
        '''
            Untere Nebendiagonale (1 Element nach unten, um 1 Elemente verkürzt)
        '''
        coeff = (damp[:,1:-1,1:-1]+rayd[1:-1,1:-1])/self.dely**2/self.param['ae']**2
        coeff = coeff + (damp[:,1:-1,1:-1]+rayd[1:-1,1:-1])*self.sn[:,1:-1]/self.cs[:,1:-1]/self.param['ae']**2/2/self.dely
        coeff = coeff - (rayd[1:-1,2:]-rayd[1:-1,:-2])/4/self.dely**2/self.param['ae']**2
        return np.pad(coeff,((0,0),(1,1),(1,1)),'constant',constant_values=(0,))
    
    
    def _a6(self,damp,rayd):
        '''
            Obere Nebendiagonale (1 Element nach rechts, um 1 Elemente verkürzt)
        '''
        coeff = (damp[:,1:-1,1:-1]+rayd[1:-1,1:-1])/self.dely**2/self.param['ae']**2
        coeff = coeff - (damp[:,1:-1,1:-1]+rayd[1:-1,1:-1])*self.sn[:,1:-1]/self.cs[:,1:-1]/self.param['ae']**2/2/self.dely
        coeff = coeff + (rayd[1:-1,2:]-rayd[1:-1,:-2])/4/self.dely**2/self.param['ae']**2
        return np.pad(coeff,((0,0),(1,1),(1,1)),'constant',constant_values=(0,))
        
        
    def _a8(self,roN2,N2,damp,alpha):
        '''
            Nebendiagonale ny Elemente nach rechts (um ny Elemente verkürzt)
        '''
        coeff = (roN2[:,1:-1,1:-1]+roN2[:,2:,1:-1])/2/self.delz**2 * (self.fcor[:,1:-1]**2 * self.ez[1:-1,:] *
                                                                        (damp[:,1:-1,1:-1]+alpha[1:-1,1:-1]))
        coeff = coeff + self.fcor[:,1:-1]**2/N2[:,1:-1,1:-1] * (alpha[2:,1:-1]-alpha[:-2,1:-1])/4/self.delz**2
        return np.pad(coeff,((0,0),(1,1),(1,1)),'constant',constant_values=(0,))
        
        
    def _a2(self,roN2,N2,damp,alpha):
        '''
            Nebendiagonale ny Elemente nach unten (um ny Elemente verkürzt)
        '''
        coeff = (roN2[:,1:-1,1:-1]+roN2[:,:-2,1:-1])/2/self.delz**2 * (self.fcor[:,1:-1]**2 * self.ez[1:-1,:] *
                                                                        (damp[:,1:-1,1:-1]+alpha[1:-1,1:-1]))
        coeff = coeff - self.fcor[:,1:-1]**2/N2[:,1:-1,1:-1] * (alpha[2:,1:-1]-alpha[:-2,1:-1])/4/self.delz**2
        return np.pad(coeff,((0,0),(1,1),(1,1)),'constant',constant_values=(0,))
        
        
    def _stack(self,data):
        '''
            stack z- and y-dimension to (z_0,y_0),(z_0,y_1), ..., (z_nz,y_ny-1),(z_nz,y_ny)
            
            data should have dimensions (nt,nz,ny)
        '''
        return np.reshape(data,(self.param['nt'],-1),order='C')
    
    def _unstack(self,data):
        return np.reshape(data,(self.param['nt'],self.param['nz'],self.param['ny']),order='C')
    
    
    def A(self,qy,roN2,N2,damp,rayd,alpha):
        '''   
            matrix A in the QG linear system
            
            CHECK THE MATRIX EQUATION
            Not sure whether I have to cut the first or the last elements
        '''
        diagonals = [self._stack(self._a2(roN2,N2,damp,alpha))[:,self.param['ny']:],
                     self._stack(self._a4(damp,rayd))[:,1:],
                     self._stack(self._a5(qy,roN2,damp,rayd,alpha))[:,:],
                     self._stack(self._a6(damp,rayd))[:,:-1],
                     self._stack(self._a8(roN2,N2,damp,alpha))[:,:-self.param['ny']],]
        offsets = [-self.param['ny'],-1,0,1,self.param['ny']]
        return diagonals, offsets
    
    
    def solve(self,A,f):
        '''
            might require a phi.todense()
        '''
        phi = []
        for i in range(self.param['nt']):
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
        f = np.zeros((self.param['nz'],self.param['ny']))
        # Meridional structure of geopotential height is set to constant 1
        # divide by coriolis to obtain streamfunction
        f[0,1:-1] = 1 / self.fcor[:,1:-1]
        return f.flatten()
    
    
    def roN2(self,N2):
        '''
            squared buoancy frequency that includes the effect of density
        '''
        result = 1 / self.ez / N2
        return result
        
        
class Diagnostic:
    def __init__(self,namelist,**kwargs):
        # define default parameters
        self.param = dict(ny=47,nz=71,nt=1,
                          top=15,bottom=0,wide=90,hem=1,
                          Ho=7000,om=7.27*10**(-5),ae=6.37*10**6,
                          s=2)
        
        # read pairs of key and value from namelist
        # ignore lines that start with '#', use whitespaces as delimited
        f = open(namelist,'r')
        for line in f:
            if line.startswith('#'):
                continue
            else:
                l = line.strip().split(' ',1)
                self.param.update({l[0]:int(l[1])})
        f.close()
        
        # update parameters from kwargs
        self.param.update(kwargs)
        
        self.delz = (self.param['top'] - self.param['bottom'])/self.param['nz']*self.param['Ho']
        width = self.param['wide'] * np.pi / 180
        self.dely = width / (self.param['ny']-1)
        
        # y-coordinate (increasing latitude)
        if self.param['hem'] > 0: # Northern Hemisphere
            y = np.linspace(0,width,self.param['ny'])
        else:
            y = np.linspace(-np.pi/2,width-np.pi/2,self.param['ny'])
        # z-coordinate [bottom,top)
        z = np.arange(self.param['bottom']*self.param['Ho'],self.param['top']*self.param['Ho'],self.delz)
        
        self.fcor = 2*self.param['om']*np.sin(y[np.newaxis,:])
        self.sn = np.sin(y[np.newaxis,:])
        self.cs = np.cos(y[np.newaxis,:])
        self.ez = np.exp(z[:,np.newaxis]/self.param['Ho'])
        
        
    def _fy(self,data):
        '''
            First y-derivative
        '''
        derv = lambda a: (a[:,:,2:] - a[:,:,:-2]) / 2 / self.dely
        result = derv(np.pad(data,((0,0),(0,0),(1,1)),mode='edge'))
        return result

    def _fyy(self,data):
        '''
            Second second y-derivative
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
            (rho/N2*f_z)_z -> expand product rule in diff
        '''
        factor = self.ez / 2 / self.delz / self.delz
        roN2 = 1 / self.ez / N2
        
        diff = lambda a,roN2: ((roN2[:,2:,:]+roN2[:,1:-1,:])*a[:,2:,:] - 
                               (roN2[:,2:,:]+2*roN2[:,1:-1,:]+roN2[:,:-2,:])*a[:,1:-1,:] + 
                               (roN2[:,1:-1,:]+roN2[:,:-2,:])*a[:,:-2,:])
        result = diff(np.pad(data,((0,0),(1,1),(0,0)),mode='edge'),np.pad(roN2,((0,0),(1,1),(0,0)),mode='edge'))
        result = result * factor
        return result
    
    
    def FN2(self,N2):
        '''
            contribution of density and N2 to the index of refraction
        '''
        f = np.sqrt(self.ez*N2)
        result = self.fcor**2 / N2 * self._fznz(f,N2)
        return result
    
    
    def l2(self,phi):
        '''
            Meridional wavenumber squared
        '''
        result = - np.real((self._fyy(phi)-self.sn/self.cs*self._fy(phi))/phi)
        result[np.isinf(result)] = 0
        result[np.isnan(result)] = 0
        return result
        
        
    def m2(self,phi,N2):
        '''
            Vertical wavenumber squared
        '''
        psizz = np.sqrt(N2/self.ez) * (self._fznz(phi,N2)-self.FN2(N2)*phi)
        result = - np.real(psizz/phi) * np.sqrt(self.ez * N2)
        result[np.isinf(result)] = 0
        result[np.isnan(result)] = 0
        return result
    
    
    def n2ref(self,phi,N2):
        '''
            Refractive index squared
        '''
        result = self.m2(phi,N2) + N2 / self.fcor**2 * self.l2(phi)
        result[np.isnan(result)] = 0
        return result
        
        
    def amplitude(self,phi):
        '''
            DO I HAVE TO DIVIDE BY FCOR to get amplitude in units of geopotential height?
        '''
        result = np.abs(phi)
        return result
    
    
    def phase(self,phi):
        result = np.arctan2(np.imag(phi),np.real(phi))
        return result
        
        
    def EPflux(self,phi,N2):
        epfy = np.imag(np.conjugate(phi)*self._fy(phi)) / self.ez / 2 / self.param['ae'] * self.param['s']
        epfz = np.imag(np.conjugate(phi)*self._fz(phi)) / self.ez / N2 * self.fcor**2  / 2 * self.param['s']
        return epfy, epfz
    
    
    def enstrophy(self,phi,N2):
        '''
            THIS IS STILL COMPLEX, where to take the real part?
        '''
        phixx = - phi * self.param['s']**2 / self.param['ae']**2 / self.cs**2
        qprime = self.fcor**2 * self._fznz(phi,N2) + phixx 
        qprime = qprime + self._fyy(phi) / self.param['ae']**2 + self._fy(phi) * self.sn / self.cs / self.param['ae']**2
        result = qprime ** 2
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