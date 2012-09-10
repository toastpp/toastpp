import numpy as np
from toast import toastmod
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
    return toastmod.ReadMesh(name)

def Write(hmesh,name):
    """Write a Toast mesh to file.

    Syntax: mesh.Write(hmesh,filename)

    Parameters:
        hmesh: mesh handle (int32)
        filename: mesh file name (string)
    """
    return toastmod.WriteMesh(hmesh,name)

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
    return toastmod.MakeMesh(nlist,elist,eltp)

def Clear(hmesh):
    """Deallocate a mesh handle.

    Syntax: mesh.Clear(hmesh)

    Parameters:
        hmesh: mesh handle

    Notes:
        After the function has returned, the mesh handle is invalid and
        should no longer be used.
    """
    return toastmod.ClearMesh(hmesh)

def Data(hmesh):
    return toastmod.MeshData(hmesh)

def SurfData(hmesh):
    return toastmod.SurfData(hmesh)

def NodeCount(hmesh):
    return toastmod.MeshNodeCount(hmesh)

def ElementCount (hmesh):
    return toastmod.MeshElementCount(hmesh)

def Dim(hmesh):
    return toastmod.MeshDim(hmesh)

def BB(hmesh):
    return toastmod.MeshBB(hmesh)

def RasterBasisPoints(hraster):
    return toastmod.RasterBasisPoints(hraster)

def RasterSolutionPoints(hraster):
    return toastmod.RasterSolutionPoints(hraster)

def ReadQM(hmesh,qmname):
    return toastmod.ReadQM(hmesh,qmname)

def ReadNim(nimname,idx=-1):
    return toastmod.ReadNim(nimname,idx)

def WriteNim(nimname,meshname,nim):
    return toastmod.WriteNim(nimname,meshname,nim)

def Sysmat(hmesh,mua,mus,ref,freq):
    rp,ci,vl=toastmod.Sysmat(hmesh,mua,mus,ref,freq)
    return sparse.csr_matrix((vl,ci,rp))

def Sysmat_CW(hmesh,mua,mus,ref,freq):
    rp,ci,vl=toastmod.Sysmat_CW(hmesh,mua,mus,ref)
    return sparse.csr_matrix((vl,ci,rp))

def Qvec(hmesh,type='Neumann',shape='Gaussian',width=1):
    rp,ci,vl=toastmod.Qvec(hmesh,type=type,shape=shape,width=width)
    m = rp.shape[0]-1
    n = toastmod.MeshNodeCount(hmesh)
    return sparse.csr_matrix((vl,ci,rp), shape=(m,n))

def Mvec(hmesh,shape='Gaussian',width=1):
    rp,ci,vl=toastmod.Mvec(hmesh,shape=shape,width=width)
    m = rp.shape[0]-1
    n = toastmod.MeshNodeCount(hmesh)
    return sparse.csr_matrix((vl,ci,rp), shape=(m,n))

def Fields(hmesh,hraster,qvec,mvec,mua,mus,ref,freq,mode='d'):
    """Calculate the nodal photon density fields.

    Syntax: phi = mesh.Fields(hmesh,raster,qvec,mvec,mua,mus,ref,freq,mode)

    Parameters:
        hmesh: mesh handle
        raster: raster handle (or -1 for mesh basis)
        qvec:   sparse matrix of source vectors
        mvec:   sparse matrix of measurement vectors
        mua:    vector of nodal absorption coefficients
        mus:    vector of nodal scattering coefficients
        ref:    vector of nodal refractive indices
        freq:   modulation frequency [MHz]
        mode:   [optional] string: 'd': direct fields only,
                'a': adjoint fields only, 'da': direct and adjoint fields

    Return values:
        phi:    Photon density fields. If mode=='da', then phi is returned as
                a tuple, where phi[0] are the direct fields, and phi[1] are
                the adjoint fields
    """
    
    if sparse.isspmatrix_csc(qvec)==True:
        qvec = sparse.csr_matrix(qvec)

    if sparse.isspmatrix_csc(mvec)==True:
        mvec = sparse.csr_matrix(mvec)
        
    return toastmod.Fields(hmesh,hraster,qvec.data,qvec.indptr,qvec.indices,mvec.data,mvec.indptr,mvec.indices,mua,mus,ref,freq,mode)


def Jacobian(hmesh,hraster,dphi,aphi,proj):
    return toastmod.Jacobian(hmesh,hraster,dphi,aphi,proj)


def Gradient(hmesh,hraster,qvec,mvec,mua,mus,ref,freq,data,sd):
    if sparse.isspmatrix_csc(qvec)==True:
        qvec = sparse.csr_matrix(qvec)

    if sparse.isspmatrix_csc(mvec)==True:
        mvec = sparse.csr_matrix(mvec)
        
    return toastmod.Gradient(hmesh,hraster,qvec.data,qvec.indptr,qvec.indices,mvec.data,mvec.indptr,mvec.indices,mua,mus,ref,freq,data,sd)


def Krylov(x,J):
    return toastmod.Krylov(x,J)


def Linesearch(x0,d,s0,p0,func):
    sl = 0
    pl = p0
    sh = s0
    x = x0 + d*sh
    ph = func(x)

    if ph < pl:
        sm = sh
        pm = ph
        sh = sh*2
        x = x0 + d*sh;
        ph = func(x)
        while ph < pm:
            sl = sm
            pl = pm
            sm = sh
            pm = ph
            sh = sh*2
            x = x0 + d*sh
            ph = func(x)

    else:
        sm = sh/2
        x = x0 + d*sm
        pm = func(x)
        while pm > pl:
            sh = sm
            ph = pm
            sm = sm/2
            x = x0 + d*sm
            pm = func(x)

    if ph < pm:
        pmin = pm
        s = sm
    else:
        a = ((pl-ph)/(sl-sh) - (pl-pm)/(sl-sm)) / (sh-sm)
        b = (pl-ph)/(sl-sh) - a*(sl+sh)
        s = -b/(2*a)
        x = x0 + d*s
        pmin = func(x)
        if pmin > pm:
            s = sm
            pmin = pm
    
    return (s,pmin)
    
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
    toastmod.Test(csrm)
