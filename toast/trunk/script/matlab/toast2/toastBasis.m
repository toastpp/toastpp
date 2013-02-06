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
            %         basis = toastBasis(..., gridfunc)
            %         basis = toastBasis(..., bb)
            %
            % Parameters:
            %         mesh [object]:
            %             mesh object to associate with the mapper
            %         bdim [integer array d]:
            %             grid dimensions
            %         gdim [integer array d]:
            %             dimensions for intermediate grid
            %         gridfunc [string]:
            %             basis functions to use for the regular grid:
            %             'LINEAR' or 'CUBIC'
            %         bb [real matrix d x 2]:
            %             bounding box for grid coverage
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
            %         By default, the regular voxel grid is sized so that
            %         it has the same size as the bounding box of the mesh.
            %         You can provide a custom bounding box with the bb
            %         parameter, e.g. to map only a subregion of the mesh.
            obj.handle = toast(uint32(26),mesh.handle,dims,varargin);
            obj.mesh = mesh;
        end
        
        function delete(obj)
            if obj.handle > 0
                toast(uint32(27),obj.handle);
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
            n = toast(uint32(71),obj.handle);
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
            b = toast(uint32(72),obj.handle);
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
            s = toast(uint32(73),obj.handle);
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
                bdim = toast(uint32(28),obj.handle);
            else
                [bdim gdim] = toast(uint32(28),obj.handle);
            end
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
            tgt = toast(uint32(29),obj.handle,mapstr,double(src));
        end
    end
    
    properties
        handle = uint64(0);
        mesh = toastMesh;
    end
end