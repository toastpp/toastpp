import numpy as np
from toast import toastmod
#import tglumpy
from scipy import sparse
from types import *

class Raster:
    """Basis mapping object.

    Syntax basis = toast.raster.Raster(mesh,grd)

    Parameters:
        mesh: toast.mesh.Mesh object containing a valid FEM mesh.
        grd:  integer array of length 2 or 3 (corresponding to mesh
              dimension) containing the grid size for the regular
              basis.
        intgrd: integer array of length 2 or 3 (corresponding to
                mesh dimension) containing the grid size for the
                intermediate basis.

    Notes:
        A raster object allows to map between an unstructured FEM
        nodal basis, and a regular grid basis. Typically this is
        used when the forward problem is defined as a FEM solver, and
        the inverse problem is solved on a regular grid.
    """

    def __init__(self, mesh, grd, intgrd=None):
        self.handle = None
        self.Make(mesh, grd, intgrd=intgrd)

    def __del__(self):
        self.Clear()
        
    def Handle(self):
        """Returns the internal raster handle.

        Syntax: handle = raster.Handle()
        """
        return self.handle

    def Make(self, mesh, grd, intgrd=None):
        """Initialise the mapper object by assigning a mesh and grid.

        Syntax: raster.Make(mesh, grd)

        Parameters:
            mesh: toast.mesh.Mesh object
            grd:  integer array of length 2 or 3 (corresponding to mesh
                  dimension) containing the grid size for the regular
                  basis.
            intgrd: integer array of length 2 or 3 (corresponding to
                    mesh dimension) containing the grid size for the
                    intermediate basis.
        """
        self.mesh = mesh
        if grd.dtype != np.int32:
            grd = np.array(grd,dtype=np.int32)
        self.grd = grd
        if intgrd is None:
            intgrd = np.copy(grd)
        elif intgrd.dtype != np.int32:
            intgrd = np.array(intgrd, dtype=np.int32)
        self.Clear()
        self.handle = toastmod.MakeRaster(mesh.handle,grd,intgrd)

    def Clear(self):
        if self.handle is not None:
            toastmod.ClearRaster(self.handle)
            self.handle = None

    def Map(self, mapstr, srcvec):
        """Map a scalar field from one basis to another.

        Syntax: tgt_coef = raster.Map(mapstr, srcvec)

        Parameters:
            mapstr: a string of the form 'S->T' defining the source
                    basis (S) to map from, and target basis (T) to map
                    to. "S" and "T" are placeholders for one of:
                    M: mesh basis
                    B: raster basis (fully populated bounding box)
                    S: raster basis (sparse; omitting voxels with no
                       mesh support)
            srcvec: array of basis coefficients in source basis

        Return values:
            tgt_coef: array of basis coefficients in target basis
        """
        if srcvec.size > max(srcvec.shape):
            raise ValueError("Only vectors are supported.")
        return toastmod.MapBasis(self.handle, mapstr, srcvec.flatten())

    def BasisPoints(self):
        return toastmod.RasterBasisPoints(self.handle)

    def SolutionPoints(self):
        return toastmod.RasterSolutionPoints(self.handle)
        
    def GLen(self):
        return toastmod.RasterGLen(self.handle)
    
    def BLen(self):
        return toastmod.RasterBLen(self.handle)
    
    def SLen(self):
        return toastmod.RasterSLen(self.handle)

    def Sol2Basis(self):
        """
        Returns a vector of length Slen() (solution dimension). Each entry contains
        the index of the corresponding element in the basis vector (range
        0 .. Blen()-1). Used to map between basis representation (full bounding
        box) and solution representation (supported|non-masked voxels only).
        """
        return toastmod.RasterSol2Basis(self.handle)

    def Basis2Sol(self):
        """
        Returns a vector of length Blen() (basis dimension). Each entry contains
        the index of the corresponding element in the solution vector (range
        0 .. Slen()-1), or -1 if the element is outside the support of the domain.
        """
        return toastmod.RasterBasis2Sol(self.handle)

    def GetMatrix(self, mapstr):
        src_key = mapstr[0]
        trg_key = mapstr[-1]
        mat_sizes = {'S': self.SLen(), 'G': self.GLen(), 'B': self.BLen(), 'M': self.mesh.NodeCount()}
        
        # Note: If basis and grid are identical, then matrices G and GI are null. Return identity mat.
        if mat_sizes[src_key] == mat_sizes[trg_key] and src_key in ['G', 'B'] and trg_key in ['G', 'B']:
            return sparse.identity(self.GLen())
        
        rp,ci,vl = toastmod.RasterMatrix(self.handle, mapstr)
        return sparse.csr_matrix((vl, ci, rp), shape=(mat_sizes[trg_key], mat_sizes[src_key]))
