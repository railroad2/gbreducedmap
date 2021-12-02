#!/usr/bin/env python3
import sys
import time

import numpy as np
import healpy as hp
import pylab as plt

import toast.tod
import toast.todmap
from astropy.io import fits
from astropy.time import Time

from gbpipe.gbparam import GBparam
from madam_pars import madam_pars
from utils import (mjd2lst, gblon, gblat, zenith_from_time, 
                   Rot_matrix_equatorial, get_boresight_vec, get_boresight_vec2, 
                   angle_from_meridian, Rotate, psi2vec_xp, vec2quat)


mpiworld, procs, rank = toast.mpi.get_world()
comm = toast.mpi.Comm(mpiworld)


def read_fits(fname):
    hdu = fits.open(fname)
    ra = hdu[1].data['ra']
    dec = hdu[1].data['dec']
    mjd = hdu[1].data['mjd']
    phase = hdu[1].data['P']
    flag = hdu[1].data['flag']
    
    ## convert data
    signal = phase2temp(phase)

    return ra, dec, mjd, signal, flag
    

def draw_hitmap(fname, nside=256):
    hdu = fits.open(fname)
    ra = hdu[1].data['ra']
    dec = hdu[1].data['dec']
    
    pix = hp.ang2pix(nside, ra, dec, lonlat=True)
    npix, nhit = np.unique(pix, return_counts=True)
    m = np.zeros(hp.nside2npix(nside))
    m[npix] = nhit

    hp.mollview(m)
    plt.show()


def detectors_from(files):
    if type(files) is str:
        files = [files]

    detectors = []
    for fn in files:
        hdu = fits.open(fn)
        obsname = hdu[0].header['obsname']
        detname = obsname.split('_')[-1]
        pix = hdu[0].header['pix']
        detectors.append(f'{detname}_{pix}')

    return detectors


def obsname_from(files):
    if type(files) is str:
        files = [files]
    
    hdu = fits.open(files[0])

    return hdu[0].header['obsname']
        

def nsample_from(files):
    if type(files) is str: 
        files = [files]

    hdu = fits.open(files[0])
    nsample = hdu[1].data.size

    return nsample


def psi_trick1(ra, dec, mjd):
    zen = zenith_from_time(mjd)
    tel = hp.ang2vec(ra, dec, lonlat=True)
    vpsi = get_boresight_vec(tel, zen)
    psi = angle_from_meridian(tel, vpsi) 

    return psi


def psi_trick2(ra, dec):
    tel = hp.ang2vec(ra, dec, lonlat=True)
    vpsi = get_boresight_vec2(tel)
    psi = angle_from_meridian(tel, vpsi) 
    
    return psi


#def test_psi_tricks(ra, dec, mjd):
def test_psi_tricks(fname):
    ra, dec, mjd, signal, flag = read_fits(fname)

    N = None
    psi1 = psi_trick1(ra[:], dec[:], mjd[:])
    psi2 = psi_trick2(ra[:], dec[:])

    print (psi1)
    print (psi2)
    print (psi1 - psi2)
    plt.plot(psi1, label='$\psi_1$ by trick 1')
    plt.plot(psi2, label='$\psi_2$ by trick 2')
    plt.plot(psi1-psi2, label='$\psi_1-\psi_2$')
    plt.xlabel('time (ms)')
    plt.ylabel('Boresight angle (rad)')
    plt.legend()

    print (np.max(ra))
    print (np.min(ra))

    zen = zenith_from_time(mjd)
    lon, lat = hp.vec2ang(zen, lonlat=True)

    print (np.max(lon))
    print (np.min(lon))

    plt.show()


def rotmat_from_angles(ra, dec, psi=None, mjd=None, deg=True):
    if psi is None:
        if mjd is None:
            psi = psi_trick2(ra, dec)
        else:
            psi = psi_trick1(ra, dec, mjd)

    rotm = Rot_matrix_equatorial(ra, dec, psi=psi, deg=deg)

    #z = (0, 0, 1)  

    #pix_sky = hp.vec2pix(nside, *(vs.T))
    #psi_sky =  

    return rotm


def pntg_detector(vdet, psidet, rotm, xp=True):
    if xp:
        pdet = psi2vec_xp([vdet], [psidet])[0]
    else:
        pdet = psi2vec([vdet], [psidet])[0]
            
    vobs = Rotate(vdet, rotm)
    pobs = Rotate(pdet, rotm)

    psis = angle_from_meridian(vobs, pobs) # in radian
    quat = np.array(vec2quat(vobs, psis))

    return quat
     

def pixpsi_detector(vdet, psidet, rotm, xp=True, nside=1024):
    if xp:
        pdet = psi2vec_xp([vdet], [psidet])[0]
    else:
        pdet = psi2vec([vdet], [psidet])[0]
            
    vobs = Rotate(vdet, rotm)
    pobs = Rotate(pdet, rotm)

    psis = angle_from_meridian(vobs, pobs) # in radian
    pixs = hp.vec2pix(nside, *(vobs.T))

    return pixs, psis


def gen_weights(psis):
    return np.array([[1]*len(psis), np.cos(2*psis), np.sin(2*psis)]).T


def define_obs(flist, pixelinfo=None): 
    if type(flist) is str:
        flist = [flist]

    obsname = obsname_from(flist)
    detectors = detectors_from(flist)
    nsamples = nsample_from(flist)
    tod = toast.tod.TODCache(comm.comm_group, detectors, nsamples)

    if pixelinfo is None:
        pixelinfo = np.genfromtxt('gbfp_20210730.txt', names=True)

    for i, fn in enumerate(flist): 
        ra, dec, mjd, signal, flag = read_fits(fn)
        dec = dec * 2

        flag = np.array(np.logical_not(flag), dtype=int)

        rotm = rotmat_from_angles(ra, dec, mjd)

        det = detectors_from(fn)[0]

        ## calculating pointings
        d, p = det.split('_')

        if d == 'GB03':
            pixelinfo_d = pixelinfo[pixelinfo['mod']==1]
            pixelinfo_dp = pixelinfo_d[pixelinfo_d['mod_pix']==int(p)]
            if len(pixelinfo_dp) == 0:
                pixelinfo_dp = pixelinfo_d[pixelinfo_d['mod_pix']==(int(p)%4+1)]
            
        elif d == 'GB04':
            pixelinfo_d = pixelinfo[pixelinfo['mod']==3]
            pixelinfo_dp = pixelinfo_d[pixelinfo_d['mod_pix']==p]

        vdet = [pixelinfo_dp['Zp_x'][0], pixelinfo_dp['Zp_y'][0], pixelinfo_dp['Zp_z'][0]]
        psidet = pixelinfo_dp['omtffr'][0]
        #pntg = pntg_detector(vdet, psidet, rotm, xp=True)
        pixs, psis = pixpsi_detector(vdet, psidet, rotm, xp=True, nside=1024)

        ## phase signal -> temperature
        temp = phase2temp(signal)

        ## put tod in TODCache
        if i == 0:
            tod.write_times(stamps=mjd)
            tod.write_common_flags(0, flags=flag)

        tod.write(det, data=temp)
        tod.write_flags(det, flags=flag)

        ## pntg
        #tod.write_pntg(det, data=pntg)

        ## pixs, pixs directly
        weights = gen_weights(psis)
        write_weights(tod, det, weights=weights)
        write_pixels(tod, det, pixels=pixs)
            

    ## define obs
    obs = {}
    obs["name"] = obsname
    obs["tod"] = tod


    return obs


def write_weights(tod, det, weights):
    nsamples = len(tod.cache['timestamps'])
    weightsname = f"weights_{det}"
    tod.cache.create(weightsname, np.float64, (nsamples, 3))
    ref = tod.cache.reference(weightsname)
    ref[:] = weights
    return


def write_pixels(tod, det, pixels):
    nsamples = len(tod.cache['timestamps'])
    pixelsname = f"pixels_{det}"
    tod.cache.create(pixelsname, np.int64, (nsamples,))
    ref = tod.cache.reference(pixelsname)
    ref[:] = pixels
    return


def phase2temp(phase):
    temp = phase * -1
    
    return temp 


def make_map(flist):
    data = toast.Data(comm)
    
    obs = define_obs(flist)
    data.obs.append(obs)

    nside = 1024
    #toast.todmap.OpPointingHpix(nside, nest=False, mode="IQU").exec(data)

    madam = toast.todmap.OpMadam()

    params = madam_pars()
    params['path_output'] = './maps/'
    params['file_root'] = 'gbtod'

    params['temperature_only'] = True
    params['write_tod'] = True
    
    madam.params = params

    madam.exec(data)

    return


def main():
    flist = sys.argv[1:]

    make_map(flist)
    
    return 

def test():
    flist = sys.argv[1]

    test_psi_tricks(flist)
    
    return 

if __name__=="__main__":
    main()
    #test()

    


