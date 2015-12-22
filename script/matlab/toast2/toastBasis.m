classdef toastBasis < handle
    % The toastBasis class represents mapping operators between the
    % unstructured finite element mesh and a regular pixel or voxel
    % grid that may be used as a reconstruction basis or for
    % visualisation purposes.
    
    methods
        function obj = toastBasis(mesh,dims,varargin)
            % Create a new basis object.
            %
            % Syntax: basis = toastBasis(mesh, bdim)
            %         basis = toastBasis(..., gdim)
            %         basis = toastBasis(..., bb)
            %         basis = toastBasis(..., gridfunc [,gridargs])
            %
            % Parameters:
            %         mesh [object]:
            %             mesh object to associate with the mapper
            %         bdim [integer array d]:
            %             grid dimensions
            %         gdim [integer array d]:
            %             dimensions for intermediate grid
            %         bb [real matrix d x 2]:
            %             bounding box for grid coverage
            %         gridfunc [string]:
            %             basis functions to use for the regular grid.
            %             Currently supported are:
            %                  LINEAR
            %                  LINEAR_V2 (experimental)
            %                  CUBIC
            %                  BLOB_GAUSS
            %                  BLOB_BESSEL
            %                  BLOB_HANNING
            %                  BLOB_RAMP
            %                  BLOB_SPLINE
            %             Note that some basis types may support additional
            %             input parameters (see notes).
            %
            % Return values:
            %         basis [object]:
            %             new basis mapper object
            %
            % Notes:  The basis mapper relies on the associated mesh to
            %         remain valid during its lifetime. Once the mesh has
            %         been deleted, the basis object should no longer be
            %         used.
            %
            %         If the gdim parameter is provided, the basis mapper
            %         creates an intermediate grid that is used for
            %         subsampling mesh values for higher accuracy.
            %         Typically, the dimensions of gdim will be a multiple
            %         of bdim, e.g. gdim = 2*bdim.
            %
            %         The shape functions used on the regular grid can be
            %         selected with the gridfunc string, to either 'LINEAR'
            %         (bi/tri-linear, default), or 'CUBIC' (bi/tri-cubic).
	    %
            %         There is also an experimental 'LINEAR_V2' map type,
            %         which is similar to 'LINEAR', but uses a more accurate
            %         mapping algorithm. It currently works only in 2D for
            %         meshes composed of 3-noded triangles. This basis type
            %         supports an additional scalar tolerance parameter for
            %         the mapping, using key 'MAPTOL'.
            %
	    %         In addition, there is a family of bases with radially
            %         symmetric basis functions (BLOB_xxx). They differ in
            %         the radial profile of the basis function. Currently
            %         supported in this group are BLOB_RAMP (linear profile)
            %         BLOB_GAUSS (truncated Gaussian), BLOB_BESSEL (Kaiser-
            %         Bessel), BLOB_HANNING (Hanning function), and
            %         BLOB_SPLINE (cubic spline profile).
	    %         All blob basis types support an additional scalar
	    %         support radius parameter, using key 'RADIUS'.
	    %         BLOB_GAUSS also supports scalar parameter 'SIGMA'.
	    %         BLOB_BESSEL also supports scalar parameter 'ALPHA'.
            %
            %         By default, the regular voxel grid is sized so that
            %         it has the same size as the bounding box of the mesh.
            %         You can provide a custom bounding box with the bb
            %         parameter, e.g. to map only a subregion of the mesh.
            obj.handle = toastmex(uint32(26),mesh.handle,dims,varargin);
            obj.mesh = mesh;
        end
        
        function delete(obj)
            if obj.handle > 0
                toastmex(uint32(27),obj.handle);
                obj.handle = 0;
                obj.mesh = toastMesh;
            end
        end
                
        function n = nlen(obj)
            % Return the dimension of the associated node basis
            %
            % Syntax: n = basis.nlen()
            %
            % Return values:
            %         n [integer]:
            %             number of nodes in the associated mesh
            %
            % See also:
            %         toastMesh.NodeCount, BLEN, SLEN
            n = toastmex(uint32(71),obj.handle);
        end
        
        function b = blen(obj)
            % Return the length of the regular grid basis
            %
            % Syntax: b = basis.blen()
            %
            % Return values:
            %         b [integer]:
            %             number of grid points
            %
            % See also:
            %         NLEN, SLEN
            b = toastmex(uint32(72),obj.handle);
        end
        
        function s = slen(obj)
            % Return the length of the masked grid basis
            %
            % Syntax: s = basis.slen()
            %
            % Return values:
            %         s [integer]:
            %             number of grid points with support overlapping
            %             the mesh volume
            %
            % Notes:  Returns the number of grid points whose area of
            %         support overlaps the problem domain defined by the
            %         mesh volume, omitting all grid points exterior to the
            %         mesh.
            %
            % See also:
            %         NLEN, BLEN
            s = toastmex(uint32(73),obj.handle);
        end
        
        function [bdim gdim] = Dims(obj)
            % Return the grid dimensions of the regular raster basis.
            %
            % Syntax: bdim = basis.Dims()
            %         [bdim gdim] = basis.Dims()
            %
            % Return values:
            %         bdim [integer array d]:
            %             grid dimensions of the regular raster basis
            %         gdim [integer array d]:
            %             grid dimensions of the intermediate basis, if
            %             applicable
            %
            % Notes:  This method returns the grid raster dimensions
            %         specified during the creation of the basis object.
            %
            % See also:
            %         TOASTBASIS
            if nargout < 2
                bdim = toastmex(uint32(28),obj.handle);
            else
                [bdim gdim] = toastmex(uint32(28),obj.handle);
            end
        end
        
	function v = Value(obj,pt,idx,idxmode)
        % Returns the value of a basis function at a given point
        %
        % Syntax: v = basis.Value(pt,idx)
        %         v = basis.Value(pt,idx,'FULLIDX')
        %         v = basis.Value(pt,coeff)
        %         v = basis.Value(pt,coeff,'NOMASK')
        %
        % Parameters:
        %         pt [real array dim]:
        %             coordinates of evaluation position
        %         idx [integer]:
        %             basis function index
        %         coeff [real array]
        %             basis coefficient vector
        %
        % Return values:
        %         v [real]:
        %             basis function value
        %
        % Notes:  By default, index 'idx' is taken to be the 'solution'
        %         basis index, i.e. excluding any basis functions with no
        %         mesh support. If the 'FULLIDX' flag is added, the index
        %         is taken to represent the 'basis' index, i.e. for the
        %         full basis grid.
        %
        %         If a basis coefficient vector is passed as the second
        %         argument, then the return value is the value of the
        %         function expressed by the basis coefficients at pt:
        %               v = sum_i b_i(pt) * coeff(i)
        %         coeff can be of length blen or slen, i.e. expressed in
        %         the B or S basis expansion.
        %         By default, points outside the mesh domain return 0. If
        %         the 'NOMASK' flag is set, the mesh mask is not applied.
        if nargin < 4
    	    v = toastmex(uint32(81),obj.handle,pt,idx);
        else
            v = toastmex(uint32(81),obj.handle,pt,idx,idxmode);
        end
    end

    function buu = Buu(obj)
        buu = toastmex(uint32(87),obj.handle);
    end
    
    function bvv = Bvv(obj)
        bvv = toastmex(uint32(88),obj.handle);
    end
    
    function buv = Buv(obj)
        buv = toastmex(uint32(89),obj.handle);
    end
    
    function duu = Duu(obj)
        duu = toastmex(uint32(98),obj.handle);
    end
    
    function dvv = Dvv(obj)
        dvv = toastmex(uint32(99),obj.handle);
    end
    
    function duv = Duv(obj)
        duv = toastmex(uint32(97),obj.handle);
    end
    
    function bvw = Bvw(obj,dim)
        bvw = toastmex(uint32(93),obj.handle,dim);
    end
    
    function sa = SupportArea(obj,idx)
        sa = toastmex(uint32(90),obj.handle,idx);
    end
    
    function Refine(obj,idx)
        toastmex(uint32(91),obj.handle,idx);
    end
    
        function tgt = Map(obj,mapstr,src)
            % Map a scalar function defined over the problem domain from
            % one basis representation to another.
            %
            % Syntax: tgt = basis.Map(mapstr,src)
            %
            % Arguments:
            %         mapstr [string]:
            %             specifies source and target basis representations
            %             (see notes)
            %         src [real or complex array]:
            %             function coefficients in source basis
            %
            % Return values:
            %         tgt [real or complex array]:
            %             function coefficients in target basis
            %
            % Notes:  The mapstr argument has the format 'X->Y', where X
            %         is a character representing the source basis, and Y a
            %         character representing the target basis. The
            %         following are supported for X and Y:
            %
            %         X/Y : basis
            %
            %         M   : mesh basis
            %         B   : regular grid basis
            %         S   : regular grid basis, excluding external points
            %         G   : intermediate grid basis
            %
            % Example:
            %         simg = basis.Map('M->S', mimg)
            %         maps function 'mimg', represented by nodal values of
            %         the unstructured mesh, to regular grid basis 'simg'
            %         excluding grid points external to the mesh.
            tgt = toastmex(uint32(29),obj.handle,mapstr,double(src));
        end
        
        function img = Sample(obj,svec,grd)
            % Samples basis coefficients on a regular grid
            %
            % Syntax: img = basis.Sample(svec,grd)
            %
            % Parameters:
            %         svec [real array slen]
            %             Vector of basis coefficients
            %         grd [integer array 2 or 3]
            %             Sampling grid sizes
            %
            % Return values:
            %         img [real array]
            %             One-dimensional array of sampled image
            %
            % Notes:  To turn the returned vector into a d-dimensional
            %         image, use img = reshape(img,grd)
            img = toastmex(uint32(83),obj.handle,svec,grd);
        end
        
        function elref = GridElref(obj)
            % Returns a list of element indices for each grid point.
            %
            % Syntax: elref = basis.GridElref()
            %
            % Return values:
            %         elref [integer array]
            %             Element index for each grid point of the basis.
            %
            % Notes:  The returned list is of length glen and contains the
            %         1-based element indices for every grid point. Grid
            %         points outside the mesh are indicated by value 0.
            elref = toastmex(uint32(63),obj.handle);
        end

	function g = ImageGradient(obj,img)
	    g = toastmex(uint32(64),obj.handle,img);
	end
    end
    
    properties
        handle = uint64(0);
        mesh = toastMesh;
    end
end
