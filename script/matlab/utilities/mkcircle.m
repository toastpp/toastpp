% Return the data for constructing a circular mesh
%
% [vtx,idx,eltp] = mkcircle(rad,nsect,nring,nbnd)
%
% rad [real]:        radius [mm]
% nsect [integer]:   number of element sectors (tangential element
%                    resolution; suggestion: 6)
% nring [integer]:   number of element rings (radial resolution, e.g. 32)
% nbnd [integer]:    number of high-resolution boundary rings (e.g. 2)
% vtx [real nx3]:    vertex list
% idx [integer ex3]: element index list
% eltp [integer e]:  element type list
%
% The returned parameters can be used directly in the toastMesh
% constructor to create a mesh instance.
%
% Example:
%   [vtx,idx,eltp] = mkcircle(25,6,32,2);
%   mesh = toastMesh(vtx,idx,eltp);

function [vtx, idx, eltp] = mkcircle (rad, nsect, nring, nbnd)
  
  eps = 1e-8;
  
  nds = 1 + (nsect * (nring-nbnd+1) * (nring-nbnd))/2 + ...
        2*nsect*(nring-nbnd)*nbnd;
  
  vtx = zeros(nds,2);

  % generate vertices
  i = 1;
  for ring=0:nring
    r = rad * ring/nring;
    if ring==0
      nrnod = 1;
    elseif ring==1
      nrnod = nsect;
    elseif nring-ring >= nbnd
      nrnod = nrnod+nsect;
    end
    for rnod=0:nrnod-1
      phi = -2*pi*rnod/nrnod;
      x = r*cos(phi);
      if abs(x) < eps
        x = 0;
      end
      y = r*sin(phi);
      if abs(y) < eps
        y = 0;
      end
      vtx(i,1) = x;
      vtx(i,2) = y;
      i = i+1;
      if nring-ring < nbnd
        phi1 = phi-pi/nrnod;
        r1 = r-0.5*rad/nring;
        x = r1*cos(phi1);
        if abs(x) < eps
          x = 0;
        end
        y = r1*sin(phi1);
        if abs(y) < eps
          y = 0;
        end
        vtx(i,1) = x;
        vtx(i,2) = y;
        i = i+1;
      end
    end
  end
  
  % generate elements
  els = nsect * (nring-nbnd) * (nring+3*nbnd);
  idx = zeros(els,3);
  
  i = 1; ri1 = 1; ro1 = 2;
  for ring=1:nring-nbnd
    ris = ri1; ros = ro1;
    rin = ro1; ron = ro1+nsect*ring;
    for sect=0:nsect-1
      idx(i,1) = ri1;
      idx(i,3) = ro1;
      if ro1+1==ron
        idx(i,2) = ros;
      else
        idx(i,2) = ro1+1;
      end
      i = i+1;
      for newel=0:ring-2
        idx(i,1) = ri1+newel;
        idx(i,3) = ro1+newel+1;
        if ri1+newel+1==rin
          idx(i,2) = ris;
        else
          idx(i,2) = ri1+newel+1;
        end
        i = i+1;
        if ri1+newel+1==rin
          idx(i,1) = ris;
        else
          idx(i,1) = ri1+newel+1;
        end
        idx(i,3) = ro1+newel+1;
        if ro1+newel+2 == ron
          idx(i,2) = ros;
        else
          idx(i,2) = ro1+newel+2;
        end
        i = i+1;
      end
      ri1 = ri1 + ring-1;
      ro1 = ro1 + ring;
    end
    if ri1==1
      ri1 = ri1+1;
    end
  end
  
  nrnod = nsect*(nring-nbnd);
  ris = ri1; ros = ro1;
  
  for n=0:nrnod-1
    idx(i,1) = ri1+n;
    idx(i,3) = ro1+2*n+1;
    if n == nrnod-1
      idx(i,2) = ris;
    else
      idx(i,2) = ri1+n+1;
    end
    i = i+1;
    idx(i,1) = ri1+n;
    idx(i,3) = ro1+2*n;
    idx(i,2) = ro1+2*n+1;
    i = i+1;
    if n == nrnod-1
      idx(i,1) = ris;
    else
      idx(i,1) = ri1+n+1;
    end
    idx(i,3) = ro1+2*n+1;
    if n == nrnod-1
      idx(i,2) = ros;
    else
      idx(i,2) = ro1+2+2*n;
    end
    i = i+1;
    idx(i,1) = ro1+2*n+1;
    idx(i,3) = ro1+2*n;
    if n==nrnod-1
      idx(i,2) = ros;
    else
      idx(i,2) = ro1+2+2*n;
    end
    i = i+1;
  end
  ri1 = ri1+nrnod;
  ro1 = ro1+nrnod*2;
  
  for ring=1:nbnd-1
    ris = ri1; ros = ro1;
    for n=0:nrnod-1
      idx(i,1) = ri1+n*2;
      idx(i,3) = ro1+2*n+1;
      if n==nrnod-1
        idx(i,2) = ris;
      else
        idx(i,2) = ri1+n*2+2;
      end 
      i = i+1;
      idx(i,1) = ri1+n*2;
      idx(i,3) = ro1+n*2;
      idx(i,2) = ro1+n*2+1;
      i = i+1;
      if n==nrnod-1
        idx(i,1) = ris;
      else
        idx(i,1) = ri1+n*2+2;
      end
      idx(i,3) = ro1+n*2+1;
      if n==nrnod-1
        idx(i,2) = ros;
      else
        idx(i,2) = ro1+n*2+2;
      end
      i = i+1;
      idx(i,1) = ro1+n*2+1;
      idx(i,3) = ro1+n*2;
      if n==nrnod-1
        idx(i,2) = ros;
      else
        idx(i,2) = ro1+n*2+2;
      end
      i = i+1;
    end
    ri1 = ri1+nrnod*2;
    ro1 = ro1+nrnod*2;
  end
  
  eltp = ones(size(idx,1),1) * 15;
