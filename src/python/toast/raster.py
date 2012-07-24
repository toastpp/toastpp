import numpy as np
from toast import toastmod
import tglumpy
from scipy import sparse
from types import *

def Make(hmesh,grd):
    return toastmod.MakeRaster(hmesh,grd)

def Clear(hraster):
    return toastmod.ClearRaster(hraster)

def MapBasis(hraster,mapstr,srcvec):
    return toastmod.MapBasis(hraster,mapstr,srcvec)

