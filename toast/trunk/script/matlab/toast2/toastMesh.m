classdef toastMesh < handle
    methods
	    function obj = toastMesh(arg1,arg2,arg3)
            if nargin > 0
                if nargin == 3
                    obj.handle = toast(uint32(1),arg1,arg2,arg3);
                else
                    obj.handle = toast(uint32(0),arg1);
                end
            end
        end

        function delete(obj)
            if obj.handle > 0
                toast(uint32(3),obj.handle);
                obj.handle = 0;
            end
        end

	    function Make(obj,vtx,idx,eltp)
            if obj.handle > 0 % Delete existing mesh
	            toast(uint32(3),obj.handle);
            end
	        obj.handle = toast(uint32(1),vtx,idx,eltp);
        end

	    function Read(obj,fname)
            if obj.handle > 0 % Delete existing mesh
	            toast(uint32(3),obj.handle);
            end
            obj.handle = toast(uint32(0),fname);
        end

	    function Write(obj,fname)
            if obj.handle > 0
	            toast(uint32(2),obj.handle,fname);
            end
        end

        function nnd = NodeCount(obj)
            if obj.handle > 0
	            nnd = toast(uint32(6),obj.handle);
            else
	            nnd = 0;
            end
	    end

        function nel = ElementCount(obj)
            if obj.handle > 0
	            nel = toast(uint32(7),obj.handle);
            else
	            nel = 0;
            end
        end

	    function dim = Dimension(obj)
            if obj.handle > 0
	            dim = toast(uint32(9),obj.handle);
            else
	            dim = 0;
            end
        end

  	    function bb = BoundingBox(obj)
	        bb = toast(uint32(8),obj.handle);
        end

        function size = Size(obj)
	        size = toast(uint32(59),obj.handle);
        end

        function [vtx,idx,eltp] = Data(obj)
	        [vtx,idx,eltp] = toast(uint32(4),obj.handle);
        end

	    function [vtx,idx,perm] = SurfData(obj)
	        [vtx,idx,perm] = toast(uint32(5),obj.handle);
        end

	    function elsize = ElementSize(obj)
	        elsize = toast(uint32(10),obj.handle);
        end

	    function [evtx,eidx,eltp] = ElementData(obj,idx)
	        [evtx,eidx,eltp] = toast(uint32(11),obj.handle,idx);
        end

	    function elid = FindElement(obj,pts)
	        elid = toast(uint32(12),obj.handle,pts);
        end

	    function elmat = Elmat(obj,idx,intstr,sideidx)
            if nargin < 4
	            elmat = toast(uint32(25),obj.handle,idx,intstr);
            else
	            elmat = toast(uint32(25),obj.handle,idx,intstr,sideidx);
            end
        end

        function ReadQM(obj,qmname)
            toast(uint32(13),obj.handle,qmname);
        end
        
	function WriteQM(obj,qmname)
	    toast(uint32(16),obj.handle,qmname);
	end

        function qvec = Qvec(obj,varargin)
            qvec = toast(uint32(18),obj.handle,varargin{:});
        end
        
        function mvec = Mvec(obj,varargin)
            mvec = toast(uint32(19),obj.handle,varargin{:});
        end
        
	    function Display(obj,varargin)
            if obj.handle > 0
                toastShowMesh(obj,varargin{:});
            end
        end
    end
    
    properties
        handle = uint64(0);
    end
end
