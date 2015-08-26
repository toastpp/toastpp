classdef toastNim < double
    methods
        function obj = toastNim(data,mname)
            if (nargin < 2 || ~ischar(mname))
                mname = '';
            end
            if nargin == 0
                val = 0;
            else
                if ischar(data)
                    [val,mname] = toastmex(uint32(22),data);
                else
                    val = double(data);
                end
            end
            obj = obj@double(val);
            obj.meshname = mname;
        end

        function Write(obj,fname,meshname)
            toastmex(uint32(56),fname,meshname,double(obj));
        end
        
        function v = Values(obj)
            v = double(obj);
        end
        
%         function b = subsref(obj,s)
%             switch s(1).type
%                 case '.'
%                     if strcmp(s(1).subs,'meshname')
%                         b = obj.meshname;
%                     end
%                 case '()'
%                     a = double(obj);
%                     b = a(s(1).subs{:});
%             end
%         end
        
%         function obj = subsasgn(obj,s,val)
%             switch s(1).type
%                 case '.'
%                     if strcmp(s.subs,'meshname')
%                         obj.meshname = val;
%                     end
%                 case '()'
%                     val = builtin('subsasgn',double(obj),s,val);
%                     obj = toastNim(val,obj.meshname);
%             end
%         end
    end
    
    properties
        meshname
    end
end
