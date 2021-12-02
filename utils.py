import numpy as np
import healpy as hp

from astropy.time import Time

from gbpipe.gbdir import (Rot_matrix_equatorial, angle_from_meridian, Rotate, 
                          angle_from_meridian, psi2vec_xp)

#GB lon, lat
gblon = -16.51028
gblat = 28.30042

def mjd2lst(mjd, lon=-16.51028):
    t = Time(mjd, format='mjd')
    lst = t.sidereal_time('mean', longitude=str(lon)+'d').degree
    return lst


def zenith_from_time(mjd):
    lst = mjd2lst(mjd)
    lat = [gblat] * len(lst)
    vec = hp.ang2vec(lst, lat, lonlat=True)
    return vec


def get_boresight_vec(telescope, zenith): 
    ## takes ~ 10 s for 86400*1000 samples
    tz = telescope - zenith 
    bs = tz.T - telescope.T * np.sum(telescope * tz, axis=1)
    norm = np.linalg.norm(bs, axis=0)
    bs = (bs / norm).T
    return bs

 
def get_boresight_vec2(telescope):
    vl = telescope[:-2]
    vr = telescope[2:]
    tg = np.array([telescope[1]-telescope[0]] + list(vr - vl) + [telescope[-1] - telescope[-2]])
    tg = (tg.T / np.linalg.norm(tg.T, axis=0)).T
    bs = np.cross(telescope, tg, axis=1)
    return bs


def vec2quat(vec, psi):
    if len(np.array(vec).shape) == 2:
        if len(np.array(psi).shape):
            quat_arr = []

            for v, p in zip(vec, psi):
                quat_arr.append(vec2quat(v, p))

            return quat_arr
        else:
            quat_arr = []

            for v in vec:
                quat_arr.append(vec2quat(v, psi))

            return quat_arr
    else:
        if vec.shape[0] != 3:
            raise Exception("Error: The length of the vector or an element of the vector array should be 3.")
            return -1;

    quat = [np.cos(psi/2)] + list(np.sin(psi/2)*vec)

    return quat


