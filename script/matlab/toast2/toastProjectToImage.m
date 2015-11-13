function p = toastProjectToImage(hProj,field)
%toastProjectToImage  - Map a volume field to a projection
%
%   Syntax: proj = toastProjectToImage(hProj, field)
%     hProj: projector handle
%     field:   field vector (NIM format)
%     proj:    projection matrix
%

  p = toastmex(uint32(2009),hProj,field);
  
