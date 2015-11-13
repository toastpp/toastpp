function p = toastProjectToField(hProj,image)
%toastProjectToField  - Map an image set to a field (transpose of image projection)
%
%   Syntax: field = toastProjectToField(hProj, image)
%     hProj: projector handle
%     field: field vector (NIM format)
%     image: image set
%

  p = toastmex(uint32(2010),hProj,image);
  
