import numpy as np
from toast import toastmod
from scipy import sparse

class Element:
    """The Element class represents a single element in a FEM mesh.
    """
    
    def __init__(self, mesh, index):
        """Constructs an element instance from a mesh and an index (>= 0)
        """
        self.mesh = mesh
        self.elid = index
        
    def Mesh(self):
        """Returns the mesh object the element is taken from.
        """
        return self.mesh

    def Index(self):
        """Returns the element's index in the mesh (>= 0)
        """
        return self.elid

    def Dof(self):
        """Returns the the array of global degree of freedom indices for the
        element nodes.
        """
        return toastmod.elementDof(self.mesh.handle, self.elid)
    
    def Size(self):
        """Returns the element size (area for 2D elements, volume for
        3D elements).
        """
        return toastmod.elementSize(self.mesh.handle, self.elid)
    
    def Region(self, val=None):
        """Returns or sets the element region index.
        """
        if val is None:
            return toastmod.elementRegion(self.mesh.handle, self.elid)
        else:
            toastmod.elementSetRegion(self.mesh.handle, self.elid, val)

        
    def Data(self):
        """Returns the geometry of the mesh element.

        Syntax evtx,eidx,eltp = element.Data()

        Return values:
            evtx: matrix of node coordinates
            eidx: list of global node indices (see element.Dof)
            eltp: element type identifier
        """
        return toastmod.elementData(self.mesh.handle, self.elid)
        
    def Mat(self, intstr, param=None, sideidx=-1):
        """Returns an integrals of a combination of shape functions and shape
        function derivatives over the surface or interior of the element.

        Syntax: mat = element.Mat(intstr)
                mat = element.Mat(intstr, param=prm)
                mat = element.Mat(intstr, sideidx=idx)

        Parameters:
            intstr:  string defining the integral (see notes)
            sideidx: side index (only for surface integrals over a single face)
            param:   nodal parameter array (only for integrals involving a
                     field parameter)

        Return value:
            mat: matrix of integral values. Each dimension of mat is of length
                 nnd, where nnd is the number of nodes in the element. The
                 number of dimensions of mat depends on the integral type. 

        Notes:
            The integral types recognised by this function are
            given by the list below, using the following naming
            convention:
            Fi: shape function for element node i,
            Di: shape function derivative for element node i,
            dFi/dxj: partial derivative of shape function for
            element node i with respect to coordinate j
            n: number of nodes associated with the element,
            d: domain dimension.
            
            int_type : integral type      : matrix dimension
            
            'F'      : \int Fi dr         :  n x 1
            'FF'     : \int Fi Fj dr      :  n x n
            'FFF'    : \int Fi Fj Fk dr   :  n x n x n
            'DD'     : \int Di Dj dr      :  n x n
            'FD'     : \int Fi Dj dr      :  n x n x d
            'FDD'    : \int Fi Dj Dk dr   :  n x n x n
            'dd'     : \int dFi/dxj dFk/dxl dr : n x d x n x d
            'BndF'   : \int Fi ds         :  n x 1
            'BndFF'  : \int Fi Fj ds      :  n x n
            
            The following integrals include an integral of a
            function defined in the element, given by nodal
            coefficients Pi, and require the 'prm' parameter to
            be specified in the call:
            
            int_type : integral type      : matrix dimension
            
            'PFF'    : \int Pi Fj Fk dr :  n x n
            'PDD'    : \int Pi Dj Dk dr :  n x n
            'BndPFF' : \int Pi Fj Fk ds :  n x n  
            
            For boundary integrals ('BndF' and 'BndFF'), the
            integrals are performed over all element sides that
            are part of the mesh surface. Alternatively, by
            specifying a second parameter (sideidx) to the call
            to Mat, the integral can be performed over a single
            side of the element (whether boundary side or not):
            
            elmat = el.Mat('BndFF',sideidx)
            
            where sideidx is the local side index (>= 1). 
            
            The returned matrix is of dimension n x n, where n is
            the number of nodes in the element. Nonzero entries
            are located at positions (i,j) where both nodes i and
            j belong to a boundary side, or to side 'sideidx', if
            applicable.
        """
        return toastmod.elementMat(self.mesh.handle, self.elid,
                                   intstr, param, sideidx)
    
