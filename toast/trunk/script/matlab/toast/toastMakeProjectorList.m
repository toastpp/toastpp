function plist = toastMakeProjectorList

%   Generates projectors (returning list of handles) for detectors specified in hMesh
%	hMesh = mesh with (external) detector positions and normals set
%	idim = Image dimensions in pixels [w h]
%	axis = Rotation axis of gantry [ax ay]
%	camType = Type of camera model = {'ORTHO', 'PINHOLE'}
%	'flen' = Focal length for PINHOLE model (mm)
%	'pixelsize' = scale for ORTHO model (mm)
