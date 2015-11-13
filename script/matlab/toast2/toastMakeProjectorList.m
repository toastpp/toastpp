function plist = toastMakeProjectorList(mesh,idim,axis,camType,varargin)

%toastMakeProjectorList - Generates projectors for detectors
%   Generates projectors (returning list of handles) for detectors specified in hMesh
%	mesh = mesh object with (external) detector positions and normals set
%	idim = Image dimensions in pixels [w h]
%	axis = Rotation axis of gantry [ax ay]
%	camType = Type of camera model = {'ORTHO', 'PINHOLE'}
%	'flen' = Focal length for PINHOLE model (mm)
%	'pixelsize' = scale for ORTHO model (mm)
%	'shift' = offset of camera position in plane of image [ix iy] in mm

  plist = toastmex(uint32(2008),mesh.handle,idim,axis,camType,varargin);
