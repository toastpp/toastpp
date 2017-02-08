classdef toastMesh < handle
    % The toastMesh class represents finite element meshes in the
    % TOAST toolbox.
    %
    % See also:
    %     toast2, toastElement
     
    methods
	    function obj = toastMesh(arg1,arg2,arg3)
            % Create a new mesh object.
            %
            % Syntax: mesh = toastMesh(vtx,idx,eltp)
            %         mesh = toastMesh(fname)
            %         mesh = toastMesh(fname,'gmsh')
            %
            % Parameters:
            %         vtx [real array n x 2 or n x 3]:
            %             list of node coordinates for n mesh nodes.
            %         idx [integer array m x s]:
            %             list of element node indices for m elements
            %         eltp [integer array m x 1]:
            %             list of element type indices (see notes)
            %         fname [string]:
            %             file name for mesh file to load
            %         'gmsh'
            %             Read from file fname in gmsh format
            %
            % Return values:
            %         mesh [object]:
            %             new toast mesh object
            %
            % Notes:  The node list is of size n x 2 for 2-D meshes and
            %         n x 3 for 3-D meshes, where n is the number of mesh
            %         nodes.
            %
            %         The element list is of size m x s, where m is the
            %         number of mesh elements, and s is the number of nodes
            %         per element. The list is 1-based. For meshes with
            %         mixed element types, where elements contain different
            %         numbers of nodes, s is the maximum number of nodes
            %         for any element. In that case, trailing index entries
            %         for elements with node numbers < s should be set to
            %         0.
            %
            %         The node order for each element in the element list
            %         must be according to the convention expected by
            %         toast. (see element header files for local node
            %         order). To test the correctness of the node order,
            %         make sure that all element sizes are reported as
            %         positive (see toastMesh/ElementSize).
            %
            %         The following element type indices are currently
            %         recognised:
            %             1: 3-noded triangle (old format, deprecated)
            %             3: 4-noded tetrahedron
            %             4: 6-noded wedge
            %             5: 8-noded regular voxel
            %             6: 6-noded triangle
            %             7: 10-noded tetrahedron
            %             8: 6-noded isoparametric triangle
            %             9: 10-noded triangle
            %            10: 10-noded isoparametric triangle
            %            11: 10-noded isoparametric tetrahedron
            %            14: 4-noded regular pixel
            %            15: 3-noded triangle
            %            16: 3-noded 3-D surface triangle
            %            17: 6-noded 3-D surface triangle
            %
            %         Boundary nodes are labelled automatically according
            %         to their mesh connectivity. To apply different
            %         boundary flags, and to add internal boundaries, use
            %         the toastMesh/MarkBoundary method.
            %
            % See also:
            %         toastMesh, ELEMENTSIZE, MAKE, READ
            if nargin > 0
                if nargin == 2 && strcmpi(arg2,'gmsh')
                    obj = toastReadGmshMesh(arg1);
                    return;
                end
                if nargin == 3
                    obj.handle = toastmex(uint32(1),arg1,arg2,arg3);
                else
                    obj.handle = toastmex(uint32(0),arg1);
                end
            end
        end

        function delete(obj)
            if obj.handle > 0
                toastmex(uint32(3),obj.handle);
                obj.handle = 0;
            end
        end

	    function Make(obj,vtx,idx,eltp)
            % Replace the geometry of a mesh.
            %
            % Syntax: mesh.Make(vtx,idx,eltp)
            %
            % Parameters:
            %         vtx [real array n x 2 or n x 3]:
            %             list of node coordinates for n mesh nodes.
            %         idx [integer array m x s]:
            %             list of element node indices for m elements
            %         eltp [integer array m x 1]:
            %             list of element type indices (see notes)
            %
            % Notes:  The current mesh geometry of the mesh object is
            %         discarded and replaced by the new mesh definition.
            %         Note that other objects which were created with a
            %         reference to the mesh object, such as basis or
            %         regularisation objects, will no longer be valid after
            %         the mesh modification.
            %
            %         The node list is of size n x 2 for 2-D meshes and
            %         n x 3 for 3-D meshes, where n is the number of mesh
            %         nodes.
            %
            %         The element list is of size m x s, where m is the
            %         number of mesh elements, and s is the number of nodes
            %         per element. The list is 1-based. For meshes with
            %         mixed element types, where elements contain different
            %         numbers of nodes, s is the maximum number of nodes
            %         for any element. In that case, trailing index entries
            %         for elements with node numbers < s should be set to
            %         0.
            %
            %         The node order for each element in the element list
            %         must be according to the convention expected by
            %         toast. (see element header files for local node
            %         order). To test the correctness of the node order,
            %         make sure that all element sizes are reported as
            %         positive (see toastMesh/ElementSize).
            %
            %         The following element type indices are currently
            %         recognised:
            %             1: 3-noded triangle (old format, deprecated)
            %             3: 4-noded tetrahedron
            %             4: 6-noded wedge
            %             5: 8-noded regular voxel
            %             6: 6-noded triangle
            %             7: 10-noded tetrahedron
            %             8: 6-noded isoparametric triangle
            %             9: 10-noded triangle
            %            10: 10-noded isoparametric triangle
            %            11: 10-noded isoparametric tetrahedron
            %            14: 4-noded regular pixel
            %            15: 3-noded triangle
            %            16: 3-noded 3-D surface triangle
            %            17: 6-noded 3-D surface triangle
            %
            %         Boundary nodes are labelled automatically according
            %         to their mesh connectivity. To apply different
            %         boundary flags, and to add internal boundaries, use
            %         the toastMesh/MarkBoundary method.
            %
            % See also:
            %         toastMesh, ELEMENTSIZE, READ
            if obj.handle > 0 % Delete existing mesh
	            toastmex(uint32(3),obj.handle);
            end
	        obj.handle = toastmex(uint32(1),vtx,idx,eltp);
        end

	    function Read(obj,fname,arg)
            % Replace the geometry of a mesh with a definition loaded from
            % a toast mesh file.
            %
            % Syntax: mesh.Read(fname)
            %         mesh.Read(fname,'gmsh')
            %
            % Parameters:
            %         fname [string]:
            %             file name for mesh file to load
            %         'gmsh'
            %             read from file fname in gmsh format
            %
            % Notes:  The current mesh geometry of the mesh object is
            %         discarded and replaced by the new mesh definition.
            %         Note that other objects which were created with a
            %         reference to the mesh object, such as basis or
            %         regularisation objects, will no longer be valid after
            %         the mesh modification.
            %
            %         Boundary nodes are labelled automatically according
            %         to their mesh connectivity. To apply different
            %         boundary flags, and to add internal boundaries, use
            %         the toastMesh/MarkBoundary method.
            %
            % See also:
            %         toastMesh
            if obj.handle > 0 % Delete existing mesh
	            toastmex(uint32(3),obj.handle);
            end
            if nargin>2
                if strcmpi(arg,'gmsh')
                    obj.handle = toastReadGmshMesh(fname);
                    return;
                end
	    end
            obj.handle = toastmex(uint32(0),fname);
        end

	    function Write(obj,fname,fmt,prm)
            % Write a toast mesh to a file.
            %
            % Syntax: mesh.Write(fname)
            %         mesh.Write(fname,fmt)
            %
            % Parameters:
            %         fname [string]:
            %             file name for mesh file to write
            %         fmt [string]:
            %             output format. Currently supported:
            %             - 'gmsh': Write mesh in Gmsh ASCII format
            %
            % Notes:  Without format string, the mesh is written in native
            %         Toast format.
            %
            %         Gmsh output currently only supports 4-noded
            %         tetrahedral meshes. Element region information is
            %         preserved as element tags, and surface meshes are
            %         also added, based on region interfaces.
            %
            % See also:
            %         toastMesh
            if obj.handle > 0
                if nargin > 2 
                    if strcmpi(fmt,'gmsh')
                        %toastmex(uint32(2),obj.handle,fname,'gmsh');
                        toastWriteGmshMesh(obj,fname);
                    elseif strcmpi(fmt,'gmsh-surf')
                        toastWriteGmshMesh(obj,fname,'surf');
                    elseif strcmpi(fmt,'gmsh-pos') && nargin > 3
                        toastWriteMeshGmshPos(obj,fname,prm);
                    end
                else
    	            toastmex(uint32(2),obj.handle,fname);
                end
            end
        end

        function WriteVtk(obj,fname,nim)
            if obj.handle > 0
                toastmex(uint32(82),obj.handle,fname,nim);
            end
        end
        
        function perm = Optimise(obj,modestr)
            % Compute a mesh node reordering.
            %
            % Syntax: perm = mesh.Optimise(mode)
            %
            % Parameters:
            %         mode [string]:
            %             MMD:     MMD reordering
            %             BANDW:   bandwidth minimisation
            %             TINNEY2: Tinney-2 reordering
            %
            % Return values:
            %         perm [integer array]:
            if obj.handle > 0
                perm = toastmex(uint32(70),obj.handle,modestr);
            end
        end
        
        function Reorder(obj,perm)
            % Re-order mesh nodes.
            %
            % Syntax: mesh.Reorder(perm)
            %
            % Parameters:
            %         perm [integer array]:
            %             permutation vector
            %
            % Notes:  On exit, the mesh node and element index lists
            %         are re-ordered according to the permutation array.
            %
            %         This can be used e.g. for mesh optimisation in
            %         combination with mesh.Optimise.
            %
            %         Nodal images created for the original mesh will need
            %         to be re-ordered as well to work with the modified
            %         mesh.
            %        
            % See also:
            %         toastMesh/Optimise
            if obj.handle > 0
                toastmex(uint32(84),obj.handle,perm);
            end
        end
        
        function el = Element(obj,elid)
            % Return an element object for one of the mesh elements.
            %
            % Syntax: el = mesh.Element(elid)
            %
            % Parameters:
            %         elid [integer]:
            %             element index (1-based)
            %
            % Return values:
            %         el [object]:
            %             toastElement object
            %
            % Notes:  The element object remains only valid as long as the
            %         mesh it refers to is valid.
            %  
            % See also:
            %         toastMesh, toastElement
            if obj.handle > 0
                el = toastElement (obj,elid);
            else
                el = 0;
            end
        end
        
        function nnd = NodeCount(obj)
            % Return the node count of a toast mesh.
            %
            % Syntax: n = mesh.NodeCount()
            %
            % Return values:
            %         n [integer]:
            %             number of nodes in the mesh
            %
            % See also:
            %         toastMesh
            if obj.handle > 0
	            nnd = toastmex(uint32(6),obj.handle);
            else
	            nnd = 0;
            end
	    end

        function nel = ElementCount(obj)
            % Return the element count of a toast mesh.
            %
            % Syntax: m = mesh.ElementCount()
            %
            % Return values:
            %         m [integer]:
            %             number of elements in the mesh
            %
            % See also:
            %         toastMesh
            if obj.handle > 0
	            nel = toastmex(uint32(7),obj.handle);
            else
	            nel = 0;
            end
        end

	    function dim = Dimension(obj)
            % Return the spatial dimension of a mesh.
            %
            % Syntax: d = mesh.Dimension()
            %
            % Return values:
            %         n [integer]:
            %             mesh dimension (2 or 3)
            %
            % Notes:  Returns 2 for a 2-D mesh, 3 for a 3-D mesh.
            %
            % See also:
            %         toastMesh
            if obj.handle > 0
	            dim = toastmex(uint32(9),obj.handle);
            else
	            dim = 0;
            end
        end

  	    function bb = BoundingBox(obj)
            % Return the mesh bounding box.
            %
            % Syntax: bb = mesh.BoundingBox()
            %
            % Return values:
            %         bb [real matrix 2 x d]:
            %             Coordinates of two opposite corners of the
            %             mesh bounding box.
            %
            % Notes:  Returns the coordinates of the lower left corner of
            %         the bounding box in bb(1,:) and of the upper right
            %         corner in bb(2,:), such that bb(1,i) <= bb(2,i),
            %         i=1..d, where d is the mesh dimension (2 or 3).
            %
            %         This function defines the bounding box by the minimum
            %         and maximum coordinates of mesh node positions.
            %         For isoparametric elements, the mesh bounding box may
            %         not be correctly represented by node positions. This
            %         is not taken into account.
            %
            % See also:
            %         toastMesh
	        bb = toastmex(uint32(8),obj.handle);
        end

        function size = Size(obj)
            % Return the size (area or volume) of a toast mesh.
            %
            % Syntax: s = mesh.Size()
            %
            % Return values:
            %         s [real]:
            %             mesh size
            %
            % Notes:  Returns the area (2-D meshes) or volume (3-D meshes)
            %         of a toast mesh. The mesh size is the sum of all
            %         element sizes.
            %
            % See also:
            %         toastMesh, ElementSize, toastElement.Size
	        size = toastmex(uint32(59),obj.handle);
        end

        function reg = Region(obj,newreg)
            % Set or retrieve the element region indices.
            %
            % Syntax: reg = mesh.Region()
            %         mesh.Region (newreg)
            %
            % Parameters:
            %         newreg [int array]:
            %             array of region indices to be assigned to the
            %             mesh elements. newreg must be of length
            %             mesh.elen()
            %
            % Return values:
            %         reg [int array]:
            %             array of current mesh element region indices
            %
            % Notes:
            %         The region index in an integer value associated with
            %         an element. It can be used to store element-specific
            %         application data, typically a region label for
            %         segmented meshes.
            %
            % See also:
            %         toastElement/Region
            if nargin > 1
                reg = toastmex(uint32(85),obj.handle,0,newreg);
            else
                reg = toastmex(uint32(85),obj.handle,0);
            end
        end
        
        function [vtx,idx,eltp] = Data(obj)
            % Return the mesh geometry.
            %
            % Syntax: [vtx,idx,eltp] = mesh.Data()
            %
            % Return values:
            %         vtx [real matrix n x d]:
            %             matrix of n d-dimensional vertices
            %         idx [integer matrix m x s]:
            %             matrix of vertex indices for each element
            %         eltp [integer array m x 1]:
            %             list of element type indices
            %
            % Notes:  Returns the node coordinates and element vertex index
            %         list for a mesh.
            %
            %         The returned vertex list is a real matrix of
            %         dimension n x d, where n is the number of nodes in
            %         the mesh, and d is the mesh dimension.
            %
            %         The returned index list is an integer matrix of
            %         dimension m x s, where m is the number of mesh
            %         elements, and s is the maximum number of nodes per
            %         element. The index list is 1-based. For elements
            %         containing fewer than s nodes, the unused entries are
            %         set to 0.
            %
            %         For a list of element type indices, see documentation
            %         of toastMesh.Make.
            %
            % See also:
            %         toastMesh, MAKE, SURFDATA
	        [vtx,idx,eltp] = toastmex(uint32(4),obj.handle);
        end

	    function [vtx,idx,perm] = SurfData(obj)
            % Return the mesh surface geometry.
            %
            % Syntax: [vtx,idx,perm] = mesh.SurfData()
            %
            % Return values:
            %         vtx [real matrix n x d]:
            %             matrix of n d-dimensional surface vertices
            %         idx [integer matrix m x s]:
            %             matrix of vertex indices for each surface element
            %         perm [integer array n x 1]:
            %             permutation array for extracting surface data
            %             from a volume array.
            %
            % Notes:  Returns the node coordinates and element vertex index
            %         list for the surface of a toast mesh.
            %
            %         The returned vertex list is a real matrix of
            %         dimension n x d, where n is the number of surface
            %         nodes, and d is the mesh dimension.
            %
            %         The returned index list is an integer matrix of
            %         dimension m x s, where m is the number of surface
            %         elements, and s is the maximum number of nodes per
            %         element.
            %
            %         The permutation array is an n-dimensional integer
            %         array that can be used to extract surface data from
            %         an array of nodal volume data, e.g.
            %         surfdata = voldata(perm)
            %
            %         The index list is 1-based. For elements containing
            %         fewer than s nodes, the unused entries are set to 0.
            %
            % See also:
            %         toastMesh, DATA
	        [vtx,idx,perm] = toastmex(uint32(5),obj.handle);
        end

        function nblist = ElementNeighbours(obj)
            % Return a list of element neighbours for each element
            %
            % Syntax: nblist = mesh.ElementNeighbours()
            %
            % Return values:
            %        nblist [integer array m x nbmax]:
            %             matrix of neighbour indices
            %
            % Note: Each row i of the returned matrix contains a list of
            %       neighbour element indices for element i.
            %
            %       An index of 0 indicates no neighbour (boundary face).
            %
            %       nbmax is the max. number of sides for any single
            %       element in the mesh. If the mesh contains mixed
            %       elements, any unused trailing entries are set to 0. 
            nblist = toastmex(uint32(100),obj.handle);
        end
        
	    function elsize = ElementSize(obj)
            % Return an array with mesh element sizes.
            %
            % Syntax: s = mesh.ElementSize()
            %
            % Return values:
            %         s [real array m x 1]:
            %             array of element sizes
            %
            % Notes:  Returns the areas (2-D meshes) or volumes (3-D
            %         meshes) of each element of a toast mesh in a column
            %         vector. Note that negative element sizes usually
            %         indicate an incorrect node order.
            %
            % See also:
            %         toastMesh, Size, toastElement.Size
	        elsize = toastmex(uint32(10),obj.handle);
        end

	    function elid = FindElement(obj,pts)
            % Find the mesh element containing a point.
            %
            % Syntax: el = mesh.FindElement(pt)
            %
            % Parameters:
            %         pt [real array d x 1]:
            %             point coordinates (2-D or 3-D coordinate array)
            %
            % Return values:
            %         el [integer]:
            %             1-based index of the element containing pt, or 0
            %             if pt is outside the mesh.
            %
            % Notes:  This function searches through the list of mesh
            %         element and checks pt against each. It can therefore
            %         be slow if applied to a large number of points.
            %
            %         For points where the element association is not
            %         unique (e.g. on element edges) this function returns
            %         the first matching element, but it is possible that
            %         no elements are found due to rounding errors.
            %
            % See also:
            %         toastMesh
	        elid = toastmex(uint32(12),obj.handle,pts);
        end

        % Return a list of outer and internal boundary faces, based on
        % element region information.
        %
        % This function extracts the faces of the outer mesh surface, and
        % any internal surfaces, and returns them as a vertex index list.
        % Internal surfaces are region boundaries, i.e. any faces that are
        % shared by elements with different region labels.
        % Outer surfaces are all faces that are not shared by two elements.
        %
        % Syntax: rsurf = mesh.RegionSurfaces
        %
        % Return values:
        %         rsurf [real nsurf x nv]:
        %             List of global vertex indices for each found surface
        %             face. Each row contains one face. nv is the max
        %             number of vertices per surface face (e.g. 3 for
        %             tetrahedral meshes).
        %
        % Notes:  This function is currently only implemented for
        %         tetrahedral meshes.
        function rsurf = RegionSurfaces(obj)
            rsurf = toastmex(uint32(94),obj.handle);
        end
            
	    function elmat = Elmat(obj,idx,intstr,sideidx)
            % Return an element matrix for a single element of a mesh.
            %
            % The returned matrix contains the integrals of products of a
            % user-specified combination of nodal shape functions and shape
            % function derivatives over the element volume or element
            % surface. The dimension of the returned matrix depends on the
            % requested integral. The element matrices are the building
            % blocks for the FEM system matrices.
            % 
            % Syntax: E = mesh.Elmat (elidx, int_type, [extra_param])
            %
            % Parameters:
            %         elidx [integer]:
            %             element index (>= 1)
            %         int_type [string]:
            %             integral type identifier (see below)
            %
            % Return values:
            %         E [real matrix]:
            %             element matrix
            %
            % Notes:  The integral types recognised by this function are
            %         given by the list below, using the following naming
            %         convention:
            %         Fi: shape function for element node i,
            %         Di: shape function derivative for element node i,
            %         dFi/dxj: partial derivative of shape function for
            %         element node i with respect to coordinate j
            %         n: number of nodes associated with the element,
            %         d: domain dimension.
            %
            %         int_type : integral type      : matrix dimension
            %
            %         'F'      : \int Fi dr         :  n x 1
            %         'FF'     : \int Fi Fj dr      :  n x n
            %         'FFF'    : \int Fi Fj Fk dr   :  n x n x n
            %         'DD'     : \int Di Dj dr      :  n x n
            %         'FD'     : \int Fi Dj dr      :  n x n x d
            %         'FDD'    : \int Fi Dj Dk dr   :  n x n x n
            %         'dd'     : \int dFi/dxj dFk/dxl dr : n x d x n x d
            %         'BndF'   : \int Fi ds         :  n x 1
            %         'BndFF'  : \int Fi Fj ds      :  n x n
            %
            %         For boundary integrals ('BndF' and 'BndFF'), the
            %         integrals are performed over all boundary sides of
            %         the element (i.e. sides forming part of the mesh
            %         surface). Alternatively, by specifying a 3th
            %         parameter (sideidx) to the call to Elmat, the
            %         integral can be performed over a single side of the
            %         element (whether boundary side or not):
            %
            %         E = mesh.Elmat(elidx,'BndFF',sideidx)
            %
            %         where sideidx is the local side index (>= 1). 
            %
            %         The returned matrix is of dimension n x n, where n is
            %         the number of nodes in the element. Nonzero entries
            %         are located at positions (i,j) where both nodes i and
            %         j belong to a boundary side, or to side 'sideidx', if
            %         applicable.
            %
            % See also:
            %         toastMesh, SYSMATCOMPONENT
            if nargin < 4
	            elmat = toastmex(uint32(25),obj.handle,idx,intstr);
            else
	            elmat = toastmex(uint32(25),obj.handle,idx,intstr,sideidx);
            end
        end

        function smat = SysmatSparsityStructure (obj)
            % Generate FEM system sparsity structure.
            %
            % Syntax: S = mesh.SysmatSparsityStructure 
            %
            % Return values:
            %         S [real sparse matrix, n x n]:
            %             empty system matrix
            %
            % Notes:
            %         Returns a sparse real n x n matrix (n: number of
            %         nodes in the mesh) with entries at locations
            %         corresponding to the non-zeros in a system matrix for
            %         this mesh. 
            %
            % See also:
            %         toastMesh, SYSMATCOMPONENT
           
            smat = toastmex(uint32(92),obj.handle);
            
        end
        
        function smat = SysmatComponent (obj,intstr,varargin)
            % Generate an FEM system matrix component.
            %
            % Syntax: S = mesh.SysmatComponent (int_tp, prm)
            %         S = mesh.SysmatComponent (int_tp, prm, 'EL')
            %
            % Parameters:
            %         int_tp [string]:
            %             integration type (see notes)
            %         prm [real column vector]:
            %             nodal parameter
            %         'EL':
            %             flag to indicate element basis
            %
            % Return values:
            %         S [complex sparse matrix, n x n]:
            %             system matrix
            %
            % Notes:
            %         Returns a sparse complex n x n matrix (n: number of
            %         nodes in the mesh) containing system matrix component
            %         S. S is constructed by integrating the parameter
            %         (prm), given by its nodal values and a combination of
            %         shape functions or shape function derivatives over
            %         all mesh elements.
            %
            %         The type of integration is selected via the int_tp
            %         string. The following are currently supported:
            %
            %         int_tp      Integral
            %         --------------------------------------------
            %         'FF'        \int Fi(r) Fj(r) dr
            %         'DD'        \int Di(r) Dj(r) dr
            %         'PFF'       \int prm(r) Fi(r) Fj(r) dr
            %         'PDD'       \int prm(r) Di(r) Dj(r) dr
            %         'BNDPFF'    \int prm(r) Fi(r) Fj(r) dr (over boundary)
            %
            %         where Fi, Di: shape function/shape function
            %         derivative for node i
            %
            %         If the 'EL' flag is present, the system matrix is
            %         calculated on an element basis, rather than a node
            %         basis. Length 'n' of all parameter vectors in that
            %         case must be equal to the number of elements.
            %         Parameters are considered constant over each element,
            %         so parameter distribution has discontinuities at
            %         element boundaries.
            %
            %         S is a building block for creating the FEM system
            %         matrix.
            %
            % Advanced usage:
            %         If a sparse matrix is passed as the last argument,
            %         this will be used as the sparsity structure for the
            %         system matrix component to be calculated. The current
            %         sparsity structure of the mesh can be extracted by
            %         calling toastMesh.SysmatSparsityStructure. Employing
            %         this technique can vastly improve performance,
            %         however caution should be excercised as no attempt is
            %         made to check the validity of the sparsity structure:
            %         improper usage will lead to undefined behaviour.
            %
            % See also:
            %         toastMesh, ELMAT, DOTSYSMAT
            
            nnd = toastmex(uint32(6),obj.handle);
            
            prm = [];
            sstruct = [];
            
            % This process is required in order not to break legacy
            % interface without optional sparsity structure.
            switch length(varargin)
                
                case 0
                    % Integral no parameter              
                case 1
                    
                    if(issparse(varargin{1}) && ismatrix(varargin{1}))
                        % Integral no parameter, with sparsity structure
                        sstruct = varargin{1};                       
                    else
                        % Integral with parameter, no sparsity strucure
                        prm = varargin{1};
                    end
                    
                case 2
                    
                    if(strcmp(varargin{2},'EL'))
                    
                        % Integral with parameter, on element, no sparsity
                        % structure
                        prm = varargin{1};
                        
                    elseif(issparse(varargin{2}) && ismatrix(varargin{2}))
                        
                        % Integral with parameter, on nodes, sparisty
                        % structure
                        prm = varargin{1};
                        sstruct = varargin{2};

                    else
                        error('mesh.SysmatComponent: unexpected arguments');
                    end
               
                case 3

                    if(issparse(varargin{3}) && ismatrix(varargin{3}))
                        % Integral with parameter, on elements, with
                        % sparsity structure
                        prm = varargin{1};
                        sstruct = varargin{3};
                    else
                        error('mesh.SysmatComponent: unexpected arguments');
                    end
                    
                otherwise
                    
                    error('mesh.SysmatComponent: unexpected arguments');
                   
            end
            
            smat = toastmex(uint32(74),obj.handle,intstr,prm,sstruct);
            
        end
        
        function mmat = Massmat (obj)
            % Generate an FEM mass matrix.
            %
            % Syntax: B = mesh.Massmat()
            %
            % Return values:
            %         B [real sparse matrix n x n]:
            %             mass matrix
            %
            % Notes:  Returns a sparse real n x n matrix (n: number of
            %         nodes in the mesh) containing mass matrix B. B is
            %         composed from element-wise integrals of products of
            %         shape functions. Unlike the system matrix S, B does
            %         not depend on any parameter distributions over the
            %         domain.
            %
            %         As an example, the mass matrix is required in the
            %         evaluation of the time-dependent FEM diffusion
            %         equation, where it appears in the time derivative
            %         term.
            %
            % See also:
            %         toastMesh, ELMAT, SYSMATCOMPONENT, dotSysmat
            mmat = toastmex(uint32(24),obj.handle);
        end
        
        function SetQM(obj,qpos,mpos)
            % Assign source and detector positions to the mesh.
            %
            % Syntax: mesh.SetQM(qpos,mpos)
            %
            % Parameters:
            %         qpos [real matrix nq x d]:
            %             array of source locations
            %         mpos [real matrix nm x d]:
            %             array of detector locations
            %
            % Notes:  This method allocates nq source and nm detector
            %         positions to the mesh. It allows the subsequent
            %         generation of right-hand side vectors with mesh.Qvec,
            %         and of measurement operators with mesh.Mvec.
            %
            %         The data link list is assumed to be dense, i.e. all
            %         sources are connected to all detectors.
            toastmex(uint32(14),obj.handle,qpos,mpos);
        end
        
        function ReadQM(obj,qmname)
            % Read a QM (source-detector) specification
            %
            % Syntax: mesh.ReadQM (fname)
            %
            % Parameters:
            %         fname [string]:
            %             QM file name
            %
            % Notes:  Reads a QM (source-detector) specification from a
            %         file into a TOAST mesh.
            %
            % See also:
            %         toastMesh, WRITEQM, DATALINKLIST, QPOS, MPOS
            toastmex(uint32(13),obj.handle,qmname);
        end
        
        function WriteQM(obj,qmname,qp,mp,dll)
            % Write source-detector definitions to a file.
            %
            % Syntax: mesh.WriteQM (qmname,qp,mp,dll)
            %
            % Parameters:
            %         qmname [string]:
            %             QM output file name
            %         qp [real array nq x dim]
            %             source positions
            %         mp [real array nm x dim]
            %             detector positions
            %         dll [sparse matrix of size nm x nq]
            %
            % Notes:  Writes the source-detector definitions associated
            %         with the mesh to a file.
            %
            %         qp contains a list of nq source coordinates,
            %         mp contains a list of nm detector coordinates.
            %         dll is a sparse matrix of size nm x nq, defining
            %         the source-detector measurement pairs. Each
            %         nonzero entry in dll indicates that the measurement
            %         for the corresponding source-detector pair should be
            %         computed.
            %
            %         If any of qp, mp or dll are set to empty [],
            %         the values stored in the mesh are used instead.
            %
            % See also:
            %         toastMesh, READQM, DATALINKLIST, QPOS, MPOS
	        toastmex(uint32(16),obj.handle,qmname,qp,mp,dll);
        end

        function linklist = DataLinkList(obj,format)
            % Returns the permutation index array for a data vector.
            %
            % Syntax: perm = mesh.DataLinkList ()
            %         perm = mesh.DataLinkList ('list')
            %         qm = mesh.DataLinkList('matrix')
            %
            % Parameters:
            %         'list':
            %             Return qm connections as a permutation array
            %             (this is the default)
            %         'matrix':
            %             Return qm connections as a sparse index matrix
            %             of dimension nm x nq
            %
            % Return values:
            %         perm [integer array]:
            %             permutation index vector
            %         qm [integer sparse matrix]:
            %             qm-connectivity matrix
            %
            % Notes:  Given a measurement setup with nq sources and nm
            %         detectors, a maximum of M = nq x nm measurements
            %         could be obtained if all detectors were recorded for
            %         every source. In practice, the measurement set is
            %         often sparse, i.e. only a subset of detectors is
            %         recorded for any given source and the total number of
            %         measurements is M' <= M.
            %
            %         DataLinkList allows to map between the full data
            %         array of size M and the sparse data array of size M'.
            %         The permutation array perm is of size M' and contains
            %         the indices in the range 1...M into the full data
            %         array.
            %
            %         Therefore
            %
            %         perm = mesh.DataLinkList();
            %         D1 = D(perm);
            %
            %         extracts the actual measurement data from a full data
            %         array D into a condensed data array D1.
            %
            %         Conversely,
            %
            %         D = zeros(n*m,1);
            %         D(perm) = D1;
            %
            %         fills the actually utilised measurements of the full
            %         data vector D with the values of the condensed vector
            %         D1, leaving the unused entries at zero.
            %
            %         If used with the 'matrix' flag, DataLinkList returns
            %         the source-detector connectivity as a sparse index
            %         matrix of size nm x nq, where a nonzero entry at
            %         qm(i,j) indicates that detector i was used with
            %         source j.
            %
            %         The nonzero entries are numbered consecutively
            %         (starting from 1) and denote the measurement position
            %         in the linear measurement vector.
            %
            %         Usage examples:
            %
            %         find(qm(:,j))
            %         Returns the detector indices used by source j. This
            %         is equivalent to row j of the QM file link list,
            %         except that indices are 1-based.
            %
            %         nonzeros(qm(:,j))
            %         Contains the indices in the linear data vector of the
            %         measurements for source j.
            %
            %         nonzeros(qm(i,:))
            %         Contains the indices in the linear data vector of the
            %         measurements for detector i.
            %
            %         data_q1 = data(nonzeros(qm(:,1)))
            %         Extracts the measurements for source 1 from linear
            %         data vector 'data'.
            %
            % See also:
            %         toastMesh, READQM, WRITEQM, QPOS, MPOS
            if nargin < 2
                format = 'list';
            elseif ~ischar(format)
                error('mesh.DataLinkList: parameter 1: string expected')
            end
            
            if strcmpi(format,'list')
    	        linklist = toastmex(uint32(17),obj.handle);
            elseif (strcmpi(format,'matrix'))
                linklist = toastmex(uint32(15),obj.handle);
            else
                error('mesh.DataLinkList: parameter 1: not recognised')
            end
        end

	    function qpos = Qpos(obj)
            % Return the source positions defined for a mesh.
            %
            % Syntax: qpos = mesh.Qpos()
            %
            % Return values:
            %         qpos [real matrix nq x d]:
            %             Array of nq source positions
            %
            % Notes:  Returns the source positions of a mesh in an nq x d
            %         matrix, where nq is the number of source positions,
            %         and d is the dimension of the mesh (d = 2 or 3).
            %         The mesh must previously have loaded a source-
            %         detector description file (see ReadQM).
            %
            % Usage example:
            %         finding the distances between sources and detectors
            %
            %         qp = mesh.Qpos();
            %         mp = mesh.Mpos();
            %         for i=1:size(qp,1)
            %             for j=1:size(mp,1)
            %                 dst(j,i) = norm(qp(i,:)-mp(j,:));
            %             end
            %         end
            %
            % See also:
            %         toastMesh, MPOS
	        qpos = toastmex(uint32(20),obj.handle);
        end

        function mpos = Mpos(obj)
            % Return the detector positions defined for a mesh.
            %
            % Syntax: mpos = mesh.Mpos()
            %
            % Return values:
            %         mpos [real matrix nm x d]:
            %             Array of nm detector positions
            %
            % Notes:  Returns the detector positions of a mesh in an nm x d
            %         matrix, where nm is the number of detector positions,
            %         and d is the dimension of the mesh (dim = 2 or 3).
            %         The mesh must previously have loaded a source-
            %         detector description file (see ReadQM).
            %
            % Usage example:
            %         finding the distances between sources and detectors
            %
            %         qp = mesh.Qpos();
            %         mp = mesh.Mpos();
            %         for i=1:size(qp,1)
            %             for j=1:size(mp,1)
            %                 dst(j,i) = norm(qp(i,:)-mp(j,:));
            %             end
            %         end
            %
            % See also:
            %         toastMesh, QPOS
            mpos = toastmex(uint32(21),obj.handle);
        end
        
        function qvec = Qvec(obj,varargin)
            % Generate a sparse matrix of source column vectors.
            %
            % Syntax: qvec = mesh.Qvec (qtype, qprof, qwidth)
            %         qvec = mesh.Qvec (prm)
            %
            % Parameters:
            %         qtype [string]:
            %             source type ('Neumann', 'Isotropic')
            %         qprof [string]:
            %             source profile ('Point', 'Gaussian', 'Cosine',
            %             'TrigBasis')
            %         qwidth [real]:
            %             source radius [mm]
            %         prm [struct]:
            %             source parameter structure (see notes)
            %
            % Return values:
            %         qvec [complex sparse matrix n x nq]:
            %             Sparse matrix containing the nodal right-hand
            %             side vector for each source location j in the
            %             j-th column.
            %
            % Notes:  Generates the RHS vectors for the FEM forward problem
            %         from source positions stored in a mesh, and source
            %         specifications passed as parameters.
            %
            %         Source specifications can be passed as a list of
            %         qtype, qprof and qwidth, or as a structure prm, which
            %         must contain the following fields:
            %         prm.type, prm.prof and prm.width.
            %  
            %         Neumann sources are defined as a diffuse boundary
            %         current (qvec is nonzero only for boundary nodes),
            %         while isotropic sources are defined in the mesh
            %         volume.
            %
            %         The width parameter defines the size of the sources
            %         (1/e radius for Gaussian sources, 1/2 radius for
            %         trigonometric sources). It is ignored for point
            %         sources.
            %
            %         The source vectors are returned as column vectors in
            %         a sparse n x nq matrix qvec, where n is the mesh node
            %         count, and nq is the number of sources. The mesh
            %         must contain source definitions (see ReadQM).
            %
            % See also:
            %         toastMesh, READQM, QPOS, MVEC
            qvec = toastmex(uint32(18),obj.handle,varargin{:});
        end
        
        function mvec = Mvec(obj,varargin)
            % Generate a sparse matrix of measurement column vectors.
            %
            % Syntax: mvec = mesh.Mvec (mprof, mwidth)
            %         mvec = mesh.Mvec (mprof, mwidth, ref)
            %
            % Parameters:
            %         mprof [string]:
            %             measurement profile ('Gaussian', 'Cosine',
            %             'TrigBasis')
            %         mwidth [real]:
            %             measurement radius [mm]
            %         ref [real, or real vector n]:
            %             refractive index used for computing the mvec
            %             operator scaling.
            %
            % Return values:
            %         mvec [complex sparse matrix n x nm]:
            %             Sparse matrix containing the nodal boundary
            %             projection operator for each detector location i
            %             in the i-th column.
            %
            % Notes:  Generates the boundary measurement operators
            %         (Direchlet-to-Neumann map) for the OT FEM forward
            %         problem from measurement specifications stored in a
            %         mesh. The vectors are returned as column vectors in a
            %         sparse n x nm matrix, where n is the mesh node count,
            %         and nm is the number of detectors. The mesh must
            %         contain measurement definitions (see ReadQM).
            %
            %         The boundary operator is normally scaled with a
            %         factor of c/2A, where c is the speed of light, and A
            %         is a term incorporating the refractive index mismatch
            %         at the boundary. This scaling factor is computed from
            %         the 'ref' term passed to Mvec (either as a vector of
            %         length n of nodal values, or as a scalar if the
            %         refractive index is homogeneous.
            %
            %         To avoid scaling with the c/2A term altogether, set
            %         the 'ref' parameter to 0.
            %
            %         Warning: if 'ref' is omitted, Mvec uses the
            %         refractive index values stored with the mesh. This
            %         behaviour may change in future versions.
            %
            % See also:
            %         toastMesh, READQM, QVEC, MPOS
            mvec = toastmex(uint32(19),obj.handle,varargin{:});
        end
        
        function intvec = IntFG(obj, f, g)
            % Integrate product of two real or complex functions over mesh.
            %
            % Syntax: mesh.IntFG(f, g)
            %
            % Parameters:
            %       f [real array n x 1]:
            %           First nodal function, real or complex double;
            %
            %       g [real array n x 1]:
            %           Second nodal function, real or complex double;
            %
            % Samuel Powell, March 2014, UCL.
            intvec = toastmex(uint32(4000),obj.handle,f,g);
        end

       
        function intvec = IntGradFGradG(obj, f, g)
            % Integrate product of the gradient of two real or complex
	    % functions over mesh.
            %
            % Syntax: mesh.IntGradFGradG(f, g)
            %
            % Parameters:
            %       f [real array n x 1]:
            %           First nodal function, real or complex double.
            %
            %       g [real array n x 1]:
            %           Second nodal function, real or complex double.
            %
            % Samuel Powell, March 2014, UCL.
            intvec = toastmex(uint32(4001),obj.handle,f,g);
        end

        function nn = NodeNeighbour (obj)
            % Generate the node neighbour matrix of the mesh.
            %
            % Syntax: nn = mesh.NodeNeighbour();
            %
            % Return values:
            %         nn [nlen x maxnbr integer matrix]
            %              Each row i of nn contains a list of indices
            %              >= 1 of neighbour nodes of i. Trailing entries
            %              of 0 are unused.
            nn = toastmex(uint32(79),obj.handle);
        end
        
        function Refine(obj,elrefine)
            % Refine elements in a mesh.
            %
            % Syntax: mesh.Refine
            %         mesh.Refine (ellist)
            %
            % Parameters:
            %         ellist [real array]
            %             list of elements to be refined
            %
            % Notes:  If no ellist is provided, all elements are refined.
            %
            %         The method also recursively refines neighbouring
            %         elements if required.
            if nargin > 1
                toastmex(uint32(77),obj.handle,elrefine);
            else
                toastmex(uint32(77),obj.handle,elrefine);
            end
        end
        
        function SplitElement(obj,elidx)
            % Split an element in the mesh.
            %
            % Syntax: mesh.SplitElement(elidx)
            %
            % Parameters:
            %         elidx
            %             element indx (>= 1)
            %
            % Notes:  currently works only with Triangle3 meshes.
            toastmex(uint32(78),obj.handle,elidx);
        end
    
        function uphase = UnwrapPhase(obj,phase,seed)
            % Unwrap a nodal phase field.
            %
            % Syntax: mesh.UnwrapPhase(phase,seed)
            %
            % Parameters:
            %         phase [real vector nlen]
            %             nodal phase field coefficients to be unwrapped
            %         seed [real vector dim]
            %             seed point for unwrapping
            %
            % Notes:  This function will fail if the target phase
            %         difference between neighbouring nodes is > pi
            uphase = toastmex(uint32(80),obj.handle,phase,seed);
        end
            
	    function h = Display(obj,varargin)
            % Display a mesh geometry and nodal field.
            %
            % Syntax: mesh.Display ()
            %         mesh.Display (field)
            %         mesh.Display (field, property, value, ...)
            %
            % Parameters:
            %         field [real array n]:
            %             nodal coefficients
            %
            % Properties/values:
            %         'range', [rmin,rmax]
            %             rescale field display to range rmin..rmax
            %             (default: determine from data range)
            %         'showgrid', true/false
            %             display the mesh grid on top of the field display
            %             (default: false)
            %         'lighting', true/false
            %             apply lighting to 3D mesh surface display
            %             (default: true; ignored for 2D mesh displays)
            %         'showcolorbar', true/false
            %             display a colorbar (default: true)
            %
            % Notes:  Displays the mesh (2D) or mesh surface (3D) and
            %         optionally a field defined by nodal coefficients
            %
            % See also:
            %         toastMesh
            if obj.handle > 0
                h = toastShowMesh(obj,varargin{:});
            end
        end
    end
    
    properties
        handle = uint64(0);
    end
end
