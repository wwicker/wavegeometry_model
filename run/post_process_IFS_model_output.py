## post_process_IFS_model_output.py
# Routines to deal with model data on hybrid coordinate model levels and irregular native grids
# - regrid to zonally averaged fields with a regular latitude and log-pressure spacing
# - resample to chosen frequency
# - use bin averaging to avoid aliasing

# wolfgang.wicker@env.ethz.ch, June 2021

import numpy as np
import xarray as xr
import os
from xhistogram.xarray import histogram as xhist
from scipy.interpolate import UnivariateSpline as Spline
from numba import guvectorize, float64


def _open_dataset(directory,chunks={},**kwargs):
    '''
        kwargs could be:
            {'combine':'nested','concat_dim':'number'}
    '''
    grb_files = [f for f in os.listdir(directory) if f.endswith('.grb')]
    grb_files.sort()
    print('\n OPEN DATASET:')
    print(grb_files)
    grb_files = [directory+f for f in grb_files]
    da = xr.open_mfdataset(grb_files,chunks=chunks,engine='cfgrib',**kwargs)
    da = da[list(da.data_vars)[0]]
    return da


def _log_pressure(lnsp,level_definition,
                 Ho=7000,ps=101300):
    '''
        pressure is in Pa
    '''
    AandB = np.loadtxt(level_definition,skiprows=1)
    hybrid = xr.DataArray(AandB[::-1,0],dims=('hybrid'),name='hybrid')
    A = xr.DataArray(AandB[::-1,1],dims=('hybrid'),coords=dict(hybrid=hybrid),name='A')
    B = xr.DataArray(AandB[::-1,2],dims=('hybrid'),coords=dict(hybrid=hybrid),name='B')
    
    pressure = A + B * np.exp(lnsp)
    z = -Ho*np.log(pressure/ps)
    return z
    

def _bin_average(data,coords,dims,bins):
    binned = xhist(*coords,bins=bins,dim=dims,weights=data)
    norm = xhist(*coords,bins=bins,dim=dims)
    binned = binned / norm
    return binned

@guvectorize(
    '(float64[:], float64[:], float64[:], float64[:])',
    '(m), (m), (n) -> (n)',
    forceobj=True
)
def _spline_interpolation(yi,xi,x,out):
    spl = Spline(xi,yi,k=5)
    out[:] = spl(x)


def _avg_time(da,axis=None,refdate='1979-01-01'):
    '''
        Average an array of np.datetime64 objects
    '''
    return ((da - np.datetime64(refdate))/10**9).mean(axis)*10**9+np.datetime64(refdate)
    
    
def _resample(ds,freq='1D'):
    try:
        ds = ds.drop(('step','time')).set_index(step='valid_time').rename(step='time')
    except ValueError:
        print('')
    new_time = ds['time'].resample(time=freq).reduce(_avg_time)
    ds = ds.resample(time=freq).mean()
    ds['time'] = new_time
    return ds


##################
def procedure1(datadir,lnspdir,level_definition,
               freq='5D',ny=321,dz=400):
    '''
    '''
    # open data variable
    da = _open_dataset(datadir,chunks=dict(step=1),combine='nested',concat_dim='number')
    
    # select Northern hemisphere
    da = da.where(da.latitude>=0,drop=True)
    
    # reindex for descending hybrid coordinate and ascending z
    da = da.reindex(hybrid=da.hybrid[::-1])
    
    # resanble in time to desired frequency
    da = _resample(da,freq)
    
    # drop last time step and upper-most level
    da = da.isel(time=slice(None,-1),hybrid=slice(None,-1))
    
    # open log of surface pressure
    lnsp = _open_dataset(lnspdir,chunks=dict(step=1),combine='nested',concat_dim='number')
    lnsp = lnsp.drop('hybrid')
    lnsp = lnsp.where(lnsp.latitude>=0,drop=True)
    lnsp = _resample(lnsp,freq)
    lnsp = lnsp.isel(time=slice(None,-1))
    
    # compute log-pressure from lnsp
    z = _log_pressure(lnsp,level_definition)
    z = z.isel(hybrid=slice(None,-1))
    
    # define vertical coordinate vector with regular gridspacing
    znew = xr.DataArray(np.arange(4900,75000+dz,dz),dims=['lev',])
    
    # 5th order spline interpolation to new z vector
    da = xr.apply_ufunc(_spline_interpolation,da,z,znew,
                        input_core_dims=[['hybrid',],['hybrid',],['lev',]],
                        output_core_dims=[['lev',]],
                        output_dtypes=[[da.dtype,]],
                        dask='parallelized')
    da = da.assign_coords(lev=znew)
    
    # compute zonal mean
    da = da.groupby('latitude').mean()
    da = da.chunk(dict(latitude=None))
    
    # define latitude vector with regular gridspacing
    lat = xr.DataArray(np.linspace(0,90,ny),dims=('lat',))
    
    # 5th order spline interpolation to new lat vector
    da = xr.apply_ufunc(_spline_interpolation,da,da.latitude,lat,
                        input_core_dims=[['latitude',],['latitude',],['lat',]],
                        output_core_dims=[['lat',]],
                        output_dtypes=[[da.dtype]],
                        dask='parallelized')
    da = da.assign_coords(lat=lat)
    
    # execute processing and load into memory
    return da.compute()