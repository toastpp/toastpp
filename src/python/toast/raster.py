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

    Notes:
        A raster object allows to map between an unstructured FEM
        nodal basis, and a regular grid basis. Typically this is
        used when the forward problem is defined as a FEM solver, and
        the inverse problem is solved on a regular grid.
    """

    def __init__(self, mesh, grd):
        self.handle = None
        self.Make(mesh, grd)

    def __del__(self):
        self.Clear()
        
    def Handle(self):
        """Returns the internal raster handle.

        Syntax: handle = raster.Handle()
        """
        return self.handle

    def Make(self, mesh, grd):
        """Initialise the mapper object by assigning a mesh and grid.

        Syntax: raster.Make(mesh, grd)

        Parameters:
            mesh: toast.mesh.Mesh object
            grd:  integer array of length 2 or 3 (corresponding to mesh
                  dimension) containing the grid size for the regular
                  basis.
        """
        self.mesh = mesh
        if grd.dtype != np.int32:
            grd = np.array(grd,dtype=np.int32)
        self.grd = grd
        self.Clear()
        self.handle = toastmod.MakeRaster(mesh.handle,grd)

    def Clear(self):
        if self.handle is not None:
            toastmod.ClearRaster(self.handle)
            self.handle = None

    def Map(self, mapstr, srcvec):
        """Map a scalar field from one basis to another.

        Syntax: tgt_coef = raster.Map(mapstr, src_coef)

        Parameters:
            mapstr: a string of the form 'S->T' defining the source
                    basis (S) to map from, and target basis (T) to map
                    to. "S" and "T" are placeholders for one of:
                    M: mesh basis
                    B: raster basis (fully populated bounding box)
                    S: raster basis (sparse; omitting voxels with no
                       mesh support)
            src_coef: array of basis coefficients in source basis

        Return values:
            tgt_coef: array of basis coefficients in target basis
        """
        return toastmod.MapBasis(self.handle, mapstr, srcvec)

    def BasisPoints(self):
        return toastmod.RasterBasisPoints(self.handle)

    def SolutionPoints(self):
        return toastmod.RasterSolutionPoints(self.handle)

