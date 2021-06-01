# Python routines to solve for quasi-geostrophic wave solutions in a given basic state
# spherical coordinate using Plumb's formulation (streamfunction=geopotential height/coriolis force)
# based on the wavegeometry model by Nili Harnik and collegues
# May 2021, Wolfgang Wicker <wolfgang.wicker@env.ethz.ch>

import numpy as np


# Define a Sponge class for all damping coefficients, define a Matrix class for the matrix equation
# all damping and matrix coefficients are computed as 3D array (time,z,y)
# stack z- and y-dimension, add all coefficients and solve the maxtrix equation A*phi=forcing
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
        return np.repeat(coeff[np.newaxis,:,:],self.param['nt'],axis=0)
    
    
    def alpha(self):
        rayd = self.rayd()

        # I don't understand the z-dependent factor
        coeff = rayd + self.param['amrad'] * (0.45*np.exp(-(self.z-50000)**2/13000**2)+0.05)[:,np.newaxis]
        return np.repeat(coeff[np.newaxis,:,:],self.param['nt'],axis=0)
        
        
    def damp(self,U):
        # get complex part of phase speed
        gamma = 1j * self.param['gammac'] * self.param['ae'] * np.cos(self.y) / self.param['s']/86400
        # get real part of phase speed from frequency which is given in days**(-1)
        cph = self.param['freq'] * self.param['ae'] * np.cos(self.y) / self.param['s']/86400
        
        coeff =  1j * self.param['s'] * (U-cph-gamma) / self.param(ae) / np.cos(self.y)
        return coeff
            
              
        
        
class Matrix:
    def __init__(self,namelist,**kwargs):
        # define default parameters
        self.param = dict(ny=47,nz=71,nt=1,
                          top=15,bottom=0,wide=90,hem=1,
                          Ho=7000,om=7.27*10**(-5),ae=6.37*10**6,
                          Ro=287,g=9.81,kappa=0.2854)
        
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
    # need to set boundary coefficients to zero,
    # np.pad(a, (1,1), 'constant', constant_values=(0, 0))
    def _a5(self,qy,roN2,damp,rayd,alpha):
        '''
            Hauptdiagonale
        '''
        
    def _a4(self,damp,rayd):
        '''
            Untere Nebendiagonale (1 Element nach unten, um 2 Elemente verkürzt)
        '''
        
    def _a6(self,damp,rayd):
        '''
            Obere Nebendiagonale (1 Element nach rechts, um 2 Elemente verkürzt)
        '''
        
        
        
    def _a8(self,roN2,N,damp,alpha):
        '''
            Nebendiagonale ny Elemente nach rechts (um 2ny Elemente verkürzt)
        '''
        
        
    def _a2(self,roN2,N,damp,alpha):
        '''
            Nebendiagonale ny Elemente nach unten (um 2ny Elemente verkürzt)
        '''
        
        
    def stack(self,data):
        # DAS IST FALSCH, NUR z- und y-axis werden gestacked
        #return data.flatten() 
       
        # stack z- and y-dimension to (z0,y0),(z0,y1), ..., (znz,yny-1),(znz,yny)
        # this 'C' order
        # np.reshape(a,(nt,-1),order='C')
        
        
        
## PROCEDURE 

# initialize model instance

# read basic state

# define sponge field similiar to the basic state fields

# stack fields

# define and solve the 2D matrix equation