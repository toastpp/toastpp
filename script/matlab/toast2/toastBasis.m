classdef toastBasis < handle
    methods
        function obj = toastBasis(hmesh,dims,varargin)
            obj.handle = toast(uint32(26),hmesh.handle,dims,varargin);
            obj.hmesh = hmesh;
        end
        
        function delete(obj)
            if obj.handle > 0
                toast(uint32(27),obj.handle);
                obj.handle = 0;
                obj.hmesh = toastMesh;
            end
        end
                
        function n = nlen(obj)
            n = toast(uint32(71),obj.handle);
        end
        
        function b = blen(obj)
            b = toast(uint32(72),obj.handle);
        end
        
        function s = slen(obj)
            s = toast(uint32(73),obj.handle);
        end
        
        function [bdim gdim] = Dims(obj)
            if nargout < 2
                bdim = toast(uint32(28),obj.handle);
            else
                [bdim gdim] = toast(uint32(28),obj.handle);
            end
        end
        
        function tgt = Map(obj,mapstr,src)
            tgt = toast(uint32(29),obj.handle,mapstr,double(src));
        end
    end
    
    properties
        handle = uint64(0);
        hmesh = toastMesh;
    end
end