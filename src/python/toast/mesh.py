import numpy as np
import ctoast
import tglumpy
from scipy import sparse
from types import *

def Read(name):
    """Read a Toast mesh from file.

    Syntax: hmesh = mesh.Read(filename)
    
    Parameters:
        filename: mesh file name (string)
        
    Return values:
        hmesh: mesh handle (int32)
    """
    return ctoast.ReadMesh(name)

def Write(hmesh,name):
    """Write a Toast mesh to file.

    Syntax: mesh.Write(hmesh,filename)

    Parameters:
        hmesh: mesh handle (int32)
        filename: mesh file name (string)
    """
    return ctoast.WriteMesh(hmesh,name)

def Make(nlist,elist,eltp):
    """Create a Toast mesh from node and element index data.

    Syntax: hmesh = mesh.Make(nlist,elist,eltp)

    Parameters:
        nlist: node coordinate list (double, n x d)
        elist: element index list (int32, m x smax)
        eltp: element type list (int32, m x 1)

    Return values:
        hmesh: mesh handle (int32)

    Notes:
        n: number of nodes, d: mesh dimension (2 or 3),
        m: number of elements, smax: max number of nodes per element
    """
    return ctoast.MakeMesh(nlist,elist,eltp)

def Clear(hmesh):
    return ctoast.ClearMesh(hmesh)

def Data(hmesh):
    return ctoast.MeshData(hmesh)

def SurfData(hmesh):
    return ctoast.SurfData(hmesh)

def NodeCount(hmesh):
    return ctoast.MeshNodeCount(hmesh)

def ElementCount (hmesh):
    return ctoast.MeshElementCount(hmesh)

def Dim(hmesh):
    return ctoast.MeshDim(hmesh)

def BB(hmesh):
    return ctoast.MeshBB(hmesh)

def MakeRaster(hmesh,grd):
    return ctoast.MakeRaster(hmesh,grd)

def ClearRaster(hraster):
    return ctoast.ClearRaster(hraster)

def RasterBasisPoints(hraster):
    return ctoast.RasterBasisPoints(hraster)

def RasterSolutionPoints(hraster):
    return ctoast.RasterSolutionPoints(hraster)

def MapBasis(hraster,mapstr,srcvec):
    return ctoast.MapBasis(hraster,mapstr,srcvec)

def ReadQM(hmesh,qmname):
    return ctoast.ReadQM(hmesh,qmname)

def ReadNim(nimname,idx=-1):
    return ctoast.ReadNim(nimname,idx)

def WriteNim(nimname,meshname,nim):
    return ctoast.WriteNim(nimname,meshname,nim)

def Sysmat(hmesh,mua,mus,ref,freq):
    rp,ci,vl=ctoast.Sysmat(hmesh,mua,mus,ref,freq)
    return sparse.csr_matrix((vl,ci,rp))

def Sysmat_CW(hmesh,mua,mus,ref,freq):
    rp,ci,vl=ctoast.Sysmat_CW(hmesh,mua,mus,ref)
    return sparse.csr_matrix((vl,ci,rp))

def Qvec(hmesh,type='Neumann',shape='Gaussian',width=1):
    rp,ci,vl=ctoast.Qvec(hmesh,type=type,shape=shape,width=width)
    return sparse.csr_matrix((vl,ci,rp))

def Mvec(hmesh,shape='Gaussian',width=1):
    rp,ci,vl=ctoast.Mvec(hmesh,shape=shape,width=width)
    return sparse.csr_matrix((vl,ci,rp))

def Fields(hmesh,hraster,qvec):
    ctoast.Fields(hmesh,hraster,qvec.data,qvec.indptr,qvec.indices)

def ShowMesh(hmesh,nim=None,col=np.array([1,1,1,1]),cmap='Grey',lighting=True,mode='Both'):
    if type(hmesh) is list:
        hm = hmesh[0]
    else:
        hm = hmesh
    if Dim(hm)==3:
        tglumpy.ShowMesh3D(hmesh,nim,col,cmap,lighting,mode)
    else:
        tglumpy.ShowMesh2D(hmesh,nim,cmap,mode)
    
def Test(csrm):
    ctoast.Test(csrm)
