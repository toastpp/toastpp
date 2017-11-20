classdef toastGmshModel < handle
methods
  function obj = toastGmshModel()
  % Create a new Gmsh model
    obj.handle = toastmex(uint32(3000));
  end

  function LoadModel(obj, fname)
    toastmex(uint32(3001), obj.handle, fname);
  end	       

  function MeshModel(obj, dim, minlength, maxlength)
      if nargin < 4
          maxlength = -1.0;
      end
      if nargin < 3
          minlength = -1.0;
      end
      toastmex(uint32(3002), obj.handle, dim, minlength, maxlength);
  end

  function RemeshModel(obj, minlength, maxlength)
    if nargin < 3
      maxlength = minlength;
    end
    toastmex(uint32(3003), obj.handle, minlength, maxlength);
  end

  function WriteMesh(obj, fname)
    toastmex(uint32(3004), obj.handle, fname);
  end
  
  function mesh = GetToastMesh(obj)
    mesh = toastMesh;
    mesh.handle = toastmex(uint32(3005), obj.handle);
  end
end

properties
  handle = uint64(0)
end

end

    
