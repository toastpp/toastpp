import numpy as np
import toast
from toast import toastmod
#import tglumpy
from scipy import sparse
from types import *

class Mesh:
    """Finite element mesh object.

    Syntax: mesh = toast.mesh.Mesh()
            mesh = toast.mesh.Mesh(name=filename)
            mesh = toast.mesh.Mesh(data=(nlist,elist,eltp))

    Parameters:
        filename: mesh file name (string)
        nlist: node coordinate list
        elist: element index list
        eltp:  element type list

    Notes:
        The default constructor creates an object without associated mesh data
        and a 0 handle.
        If name is provided, the mesh is initialised from file.
        If the data tuple is provided, the mesh is initialised from that.
        See mesh.Make() for details of the nlist, elist, eltp parameters.
        name and data should not both be assigned.
    """
    
    def __init__(self,name=None,data=None):
        self.handle = 0
        if name is not None:
            self.Read(name)
        elif data is not None:
            self.Make(data[0],data[1],data[2])

    def __del__(self):
        self.Clear()
        
    def Handle(self):
        """Returns the internal mesh handle.
        
        Syntax: handle = mesh.Handle()
        """
        return self.handle

    def Read(self, name):
        """Read a Toast mesh from file.

        Syntax: handle = mesh.Read(name)
    
        Parameters:
            name: mesh file name (string)
        
        Return values:
            handle: mesh handle (int32). Usually, the returned handle will not
            be required.
        """
        self.handle = toastmod.ReadMesh(name)
        return self.handle

    def Write(self, name):
        """Write a Toast mesh to file.

        Syntax: mesh.Write(name)

        Parameters:
            name: mesh file name (string)
        """
        toastmod.WriteMesh(self.handle, name)

    def Make(self, nlist, elist, eltp):
        """Create a Toast mesh from node and element index data.
        
        Syntax: handle = mesh.Make(nlist,elist,eltp)

        Parameters:
            nlist: node coordinate list (double, n x d)
            elist: element index list (int32, m x smax)
            eltp: element type list (int32, m x 1)

        Return values:
            handle: mesh handle (int32). Usually the returned handle will not
            be required.

        Notes:
            For eltp values, see <toastroot>/src/libfe/element.h:ELID_*
            n: number of nodes, d: mesh dimension (2 or 3),
            m: number of elements, smax: max number of nodes per element
        """
        self.handle = toastmod.MakeMesh(nlist, elist, eltp)

    def Clear(self):
        """Deallocate the associated mesh data
        
        Syntax: mesh.Clear()

        Notes:
            Releases the mesh data and resets the handle to 0
        """
        if self.handle != 0:
            toastmod.ClearMesh(self.handle)
            self.handle = 0

    def Data(self):
        """Returns node, element and element type data for the mesh

        Syntax: nlist,elist,eltp = mesh.Data()

        Return values:
            nlist: node coordinate list (double, n x d)
            elist: element index list (int32, m x smax)
            eltp: element type list (int32, m x 1)

        Notes:
            If no mesh has been loaded or constructed for the object, the
            return value is None.
        """
        return toastmod.MeshData(self.handle)

    def SurfData(self):
        """Returns node and element lists for the mesh surface faces.

        Syntax: nlist,elist,perm = mesh.SurfData()

        Return values:
            nlist: node coordinate list (double, n x d)
            elist: element index list (int32, m x smax)
            perm: permutation array (int32, n x 1)

        Notes:
            Returns the node coordinates and element vertex index list for the
            mesh surface.

            The returned vertex list is a real matrix of dimension n x d, where
            n is the number of surface nodes, and d is the mesh dimension.

            The returned index list is an integer matrix of dimension m x smax,
            where m is the number of surface elements, and smax is the maximum
            number of nodes per element.

            The permutation array is an n-dimensional integer array that can be
            used to extract surface data from an array of nodal volume data,
            e.g. surfdata = voldata[perm]
        """
        return toastmod.SurfData(self.handle)

    def NodeCount(self):
        """Returns the number of mesh nodes.
        """
        return toastmod.MeshNodeCount(self.handle)
    
    def ElementCount (self):
        """Returns the number of mesh elements.
        """
        return toastmod.MeshElementCount(self.handle)

    def Dim(self):
        """Returns the mesh dimension (2 or 3).
        """
        return toastmod.MeshDim(self.handle)

    def BoundingBox(self):
        """Returns the extents of the bounding box of the mesh.
        """
        return toastmod.MeshBB(self.handle)

    def ElementSize(self):
        """Returns an array of element sizes.
        """
        return toastmod.elementSize(self.handle, -1)

    def ElementRegion(self, val=None):
        """Returns or sets the region indices of all mesh elements.

        Syntax: reg = mesh.ElementRegion()
                mesh.ElementRegion(reg)

        reg: int array of length elen (number of elements) containing
        the region indices.
        """
        if val is None:
            return toastmod.elementRegion(self.handle, -1)
        else:
            toastmod.meshSetRegion(self.handle, val)
            
    def Element(self, elid):
        """Returns the element for element index elid (>= 0)
        """
        return toast.Element(self, elid)
        
    def ReadQM(self, name):
        return toastmod.ReadQM(self.handle, name)

    def Qvec(self, type='Neumann', shape='Gaussian', width=1):
        rp,ci,vl=toastmod.Qvec(self.handle, type=type, shape=shape, width=width)
        n = rp.shape[0]-1
        m = self.NodeCount()
        return sparse.csc_matrix((vl,ci,rp), shape=(m,n))

    def Mvec(self, shape='Gaussian', width=1, ref=1):
        rp,ci,vl=toastmod.Mvec(self.handle, shape=shape, width=width, ref=ref)
        n = rp.shape[0]-1
        m = self.NodeCount()
        return sparse.csc_matrix((vl,ci,rp), shape=(m,n))

    def ReadNim(self, name, idx=-1):
        nim = toastmod.ReadNim(name, idx)
        assert nim.size == self.NodeCount(), "Nim length doesn't match node count"
        return nim
    
    def Sysmat(self, mua, mus, ref, freq):
        """Returns the FEM stiffness matrix for the DOT problem.

        Parameters:
            mua: nodal absorption values [1/mm]
            mus: nodal scattering values [1/mm]
            ref: nodal refractive index values
            freq: modulation frequency [MHz], (0 for CW problem)

        Return value:
            Complex sparse stiffness matrix in CSR format.
        """
        rp,ci,vl=toastmod.Sysmat(self.handle, mua, mus, ref, freq)
        return sparse.csr_matrix((vl,ci,rp))

    def Fields(self, raster, qvec, mua, mus, ref, freq):
        """Calculate the nodal photon density fields.
        
        Syntax: phi = mesh.Fields(raster,qvec,mua,mus,ref,freq)

        Parameters:
            raster: raster object (or None for mesh basis)
            qvec:   sparse matrix of source vectors
            mua:    vector of nodal absorption coefficients
            mus:    vector of nodal scattering coefficients
            ref:    vector of nodal refractive indices
            freq:   modulation frequency [MHz]

        Return values:
            phi:    Matrix of photon density fields, one column per source.

        Note:
            If a valid raster object is passed, the resulting fields will
            automatically be mapped to that basis. Otherwise, the fields are
            returned in the mesh basis.
        """
        if sparse.isspmatrix_csr(qvec)==True:
            qvec = sparse.csc_matrix(qvec)
            # convert to per-source order

        nq = qvec.shape[1]

        if raster is None:
            rhandle = -1
        else:
            rhandle = raster.handle
            
        return toastmod.Fields(self.handle, rhandle,
            nq, qvec.data, qvec.indptr, qvec.indices,
            mua, mus, ref, freq)

    def Jacobian(self, hraster, dphi, aphi, proj):
        """Generate the Jacobian matrix for the frequency domain DOT problem.

        Syntax: J = mesh.Jacobian(raster,dphi,aphi,proj)

        Parameters:
            raster: raster handle (or -1 for mesh basis)
            dphi: complex matrix of nodal direct fields (one column per source)
            aphi: complex matrix of nodal adjoint fields (one column per detector)
            proj: complex vector of projections (forward data)

        Return value:
            J: Jacobian matrix (real, dense)

        Note:
            Calculates the derivative of the data (log amplitude and phase) with
            respect to the coefficients (absorption and diffusion) of the
            forward operator.

            J consists of 4 blocks: d lnmod / d mua (top left), d lnmod / d kappa
            (top right), d phase / d mua (bottom left), and d phase / d kappa
            (bottom right). Dimension of J is 2m x 2n, where m is the number of
            measurements, and n is the dimension of the inverse basis.

            If raster is set to -1, the Jacobian is constructed directly in the
            mesh basis. n in that case is equal to the number of nodes.
        """
        return toastmod.Jacobian(self.handle, hraster, dphi, aphi, proj)

    def Show(self, nim=None, col=np.array([1,1,1,1]), cmap='Grey',
             lighting=True, mode='Both'):
        if self.Dim() == 3:
            tglumpy.ShowMesh3D(self.handle, nim, col, cmap, lighting, mode)
        else:
            tglumpy.ShowMesh2D(self.handle, nim, cmap, mode)
    

        
def ReadNim(nimname,idx=-1):
    return toastmod.ReadNim(nimname,idx)

def WriteNim(nimname,meshname,nim):
    return toastmod.WriteNim(nimname,meshname,nim)

def Sysmat_CW(hmesh,mua,mus,ref,freq):
    rp,ci,vl=toastmod.Sysmat_CW(hmesh,mua,mus,ref)
    return sparse.csr_matrix((vl,ci,rp))

def Jacobian(hmesh,hraster,dphi,aphi,proj):
    return toastmod.Jacobian(hmesh,hraster,dphi,aphi,proj)


def Gradient(hmesh,hraster,qvec,mvec,mua,mus,ref,freq,data,sd):
    if sparse.isspmatrix_csr(qvec)==True:
        qvec = sparse.csc_matrix(qvec)

    if sparse.isspmatrix_csr(mvec)==True:
        mvec = sparse.csc_matrix(mvec)
        
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
    
def Test(csrm):
    toastmod.Test(csrm)
