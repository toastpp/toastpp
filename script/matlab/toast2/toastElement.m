classdef toastElement < handle
    % The toastElement class represents an element inside a finite element
    % mesh in the TOAST toolbox.
    %
    % Element objects remain valid only as long as the mesh object they
    % refer to remains valid.
    %
    % See also:
    %     toast2, toastMesh
    
    methods
        function obj = toastElement(mesh,ind)
            % Create a new element object.
            %
            % Syntax: el = toastElement(mesh,elid)
            %
            % Parameters:
            %         mesh [object]:
            %             mesh object
            %         elid [integer]:
            %             element index (1-based)
            %
            % Return values:
            %         el [object]:
            %             new element object
            %
            % See also:
            %         toastElement
            obj.mesh = mesh;
            obj.elid = ind;
        end
        
        function dof = Dof(obj)
            % Return the global degrees of freedom of the element nodes.
            %
            % Syntax: dof = el.Dof()
            %
            % Return values:
            %         dof [integer vector]:
            %             global degree of freedom for each element node.
            %
            % See also:
            %         toastElement
            dof = toastmex(uint32(75),obj.mesh.handle,obj.elid);
        end
        
        function size = Size(obj)
            % Return the element size.
            %
            % Syntax: s = el.Size()
            %
            % Return values:
            %         s [real]:
            %             element size
            %
            % Notes:  Returns the area (2-D meshes) or volume (3-D meshes)
            %         of the element. Note that a negative size usually
            %         indicates an incorrect node order.
            % See also:
            %         toastElement, toastMesh.ElementSize
            size = toastmex(uint32(10),obj.mesh.handle,obj.elid);
        end
        
        function reg = Region(obj,newreg)
            % Set or retrieve the element region index.
            %
            % Syntax: reg = el.Region()
            %         el.Region (newreg)
            %
            % Parameters:
            %         newreg [int]:
            %             new region index to be assigned to the element
            %
            % Return values:
            %         reg [int]:
            %             current element region index
            %
            % Notes:
            %         The region index in an integer value associated with
            %         an element. It can be used to store element-specific
            %         application data, typically a region label for
            %         segmented meshes.
            %
            % See also:
            %         toastMesh/Region
            if nargin > 1
                reg = toastmex(uint32(85),obj.mesh.handle,obj.elid,newreg);
            else
                reg = toastmex(uint32(85),obj.mesh.handle,obj.elid);
            end
        end
        
	    function [evtx,eidx,eltp] = Data(obj)
            % Return the geometry of the mesh element.
            %
            % Syntax: [evtx,eidx,eltp] = el.Data()
            %
            % Return values:
            %         evtx [real matrix en x d]:
            %            list of element node coordinates
            %         eidx [integer matrix em x s]:
            %            list of global node indices
            %         eltp [integer]:
            %            element type index
            %
            % Notes:  Returns the node coordinates, global node indices
            %         and element type for the element.
            %
            %         For a list of element type indices, see
            %         toastMesh/Make.
            %
            % See also:
            %         toastElement, toastMesh.Make
	        [evtx,eidx,eltp] = toastmex(uint32(11),obj.mesh.handle,obj.elid);
        end

        function elmat = Mat(obj,intstr,sideidx)
            % Return an element matrix for a single element of a mesh.
            %
            % The returned matrix contains the integrals of products of a
            % user-specified combination of nodal shape functions and shape
            % function derivatives over the element volume or element
            % surface. The dimension of the returned matrix depends on the
            % requested integral. The element matrices are the building
            % blocks for the FEM system matrices.
            % 
            % Syntax: elmat = el.Mat (int_type)
            %         elmat = el.Mat (ind_type, prm)
            %         elmat = el.Mat (ind_type, sideidx)
            %         elmat = el.Mat (ind_type, prm, sideidx)
            %
            % Parameters:
            %         int_type [string]:
            %             integral type identifier (see below)
            %         prm [real array]:
            %             nodal parameter values on the element nodes
            %         sideidx [integer]:
            %             local side index (for boundary integrals)
            %
            % Return values:
            %         elmat [real matrix]:
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
            %         'Fdd'    : \int Fi dFj/dxk dFl/dxm dr :  n x n x d x n x d
            %         'BndF'   : \int Fi ds         :  n x 1
            %         'BndFF'  : \int Fi Fj ds      :  n x n
            %
            %         The following integrals include an integral of a
            %         function defined in the element, given by nodal
            %         coefficients Pi, and require the 'prm' parameter to
            %         be specified in the call:
            %
            %         int_type : integral type      : matrix dimension
            %
            %         'PFF'    : \int Pi Fj Fk dr :  n x n
            %         'PDD'    : \int Pi Dj Dk dr :  n x n
            %         'Pd'     : \int Pi dFj/dxk dr : n x d
            %         'Pdd'    : \int Pi dFj/dxk dFl/dxm dr : n x d x n x d
            %         'BndPFF' : \int Pi Fj Fk ds :  n x n  
            %
            %         For boundary integrals ('BndF' and 'BndFF'), the
            %         integrals are performed over all element sides that
            %         are part of the mesh surface. Alternatively, by
            %         specifying a second parameter (sideidx) to the call
            %         to Mat, the integral can be performed over a single
            %         side of the element (whether boundary side or not):
            %
            %         elmat = el.Mat('BndFF',sideidx)
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
            %         toastElement, toastMesh.SysmatComponent
            if nargin < 3
                elmat = toastmex(uint32(25),obj.mesh.handle,obj.elid,intstr);
            else
                elmat = toastmex(uint32(25),obj.mesh.handle,obj.elid,intstr,sideidx);
            end
        end
    
        function fun = ShapeFun(obj,pt,frame)
            % Return element shape functions at multiple points
            %
            % Syntax: fun = el.ShapeFun(pts)
            %         fun = el.ShapeFun(pts,frame)
            %
            % Parameters:
            %         pts [real matrix d x np]:
            %             array of evaluation point coordinates
            %         frame [string]:
            %             reference frame for points pts: 'local' (default)
            %             or 'global'.
            %
            % Return values:
            %         fun [real matrix nn x np]:
            %             array of element shape functions
            %
            % Notes:  For each of the supplied evalulation points, this
            %         method returns the values of the shape functions
            %         associated with all element nodes.
            %
            %         Parameter 'pts' is a matrix of dimension d x np,
            %         representing list of np evaluation points pts(:,i).
            %         d is the problem dimension (2 or 3).
            %
            %         By default, input point coordinates are interpreted
            %         in the local element frame. If passing coordinates in
            %         the global mesh frame, use flag 'global'.
            %
            %         The dimension of the returned matrix is nn x np,
            %         where np is the number of input points, and nn is the
            %         number of element nodes. Thus, fun(i,j) is the
            %         shape function for the j-th point and i-th node. 
            %
            %         The order of the shape function array is defined by
            %         the local element node order.
            %
            % See also:
            %         toastElement, SHAPEDER
            if nargin < 3
                fun = toastmex(uint32(57),obj.mesh.handle,obj.elid,pt);
            else
                fun = toastmex(uint32(57),obj.mesh.handle,obj.elid,pt,frame);
            end
        end

        function der = ShapeDer(obj,pt,frame)
            % Return element shape function derivatives at multiple points
            %
            % Syntax: der = el.ShapeDer(pts)
            %         der = el.ShapeDer(pts,frame)
            %
            % Parameters:
            %         pt [real matrix np x d]:
            %             array of evaluation point coordinates
            %         frame [string]:
            %             reference frame for points pts: 'local' (default)
            %             or 'global'.
            %
            % Return values:
            %         der [real vector nn x d x np]:
            %             array of element shape function derivatives
            %
            % Notes:  For each of the supplied evaluation points, this
            %         method returns the values of the shape function
            %         derivatives associated with all element nodes.
            %
            %         Parameter 'pts' is a matrix of dimension np x d,
            %         representing list of np evaluation points pts(i,:).
            %         d is the problem dimension (2 or 3).
            %
            %         By default, input point coordinates are interpreted
            %         in the local element frame. If passing coordinates in
            %         the global mesh frame, use flag 'global'.
            %
            %         The dimension of the returned matrix is nn x d x np,
            %         where np is the number of input points, d is the
            %         problem dimension (2 or 3), and nn is the number of
            %         element nodes. Thus, der(i,j,k) is the shape
            %         function derivative for the k-th point, i-th node in
            %         the j-th axis direction.
            %
            %         The order of the shape function array is defined by
            %         the local element node order.
            %
            % See also:
            %         toastElement, SHAPEFUNC
            if nargin < 3
                der = toastmex(uint32(58),obj.mesh.handle,obj.elid,pt);
            else
                der = toastmex(uint32(58),obj.mesh.handle,obj.elid,pt,frame);
            end
        end
    end
    
    properties
        % Reference to a toastMesh object.
        % See also: toastMesh
        mesh = 0;
        
        % Element index (>= 1)
        elid = 0;
    end
end
