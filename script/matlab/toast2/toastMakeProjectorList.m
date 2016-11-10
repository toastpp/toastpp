function plist = toastMakeProjectorList(mesh,idim,axis,camType,varargin)

%toastMakeProjectorList - Generates projectors for detectors
%
%   Generates projectors (returning list of handles) for detectors specified in mesh
%
%   toastMakeProjectorList(mesh,idim,axis,camType)
%   mesh:      mesh object with (external) detector positions and normals set
%   idim:      Image dimensions in pixels [w h]
%   axis:      Rotation axis of gantry [ax ay]
%   camType:   Type of camera model = {'ORTHO', 'PINHOLE'}
%
%   toastMakeProjectorList(..., 'pixelsize',pixelsize)
%   pixelsize: pixel width and height (mm)
%
%   toastMakeProjectorList(..., 'flen',flen)
%   flen:      Focal length for PINHOLE model (mm)
%
%   toastMakeProjectorList(..., 'shift',shift)
%   shift:     offset of camera position in plane of image [ix iy] in mm
%
%   toastMakeProjectorList(..., 'driver',driver)
%   driver:    'software' (default) or 'mesagl'. The MesaGL driver is only
%              available if Toast was compiled with Mesa-3D support.

  plist = toastmex(uint32(2008),mesh.handle,idim,axis,camType,varargin);
