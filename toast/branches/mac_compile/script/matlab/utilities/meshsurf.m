function f=meshsurf (hMesh,prm)
  [vtx,idx,perm] = toastSurfData(hMesh);
  [pmin pmax] = toastMeshBB(hMesh);
  msize=pmax-pmin;
  scale = norm(msize)/100;
  [sx,sy,sz] = sphere(10);
  sx = sx.*scale;
  sy = sy.*scale;
  sz = sz.*scale;
  [sf, sv, sc] = surf2patch(sx,sy,sz);
  sc = zeros(length(sv),3); sc(:,1)=1;
  nvtx = size(vtx,1);
  col = zeros(nvtx,3); col(:,2)=1;
  if nargin > 1
      sprm = prm(perm);
      if max(sprm) > min(sprm)
          col(:,2) = (sprm-min(sprm))/(max(sprm)-min(sprm));
      end
  end
  qp=toastQPos(hMesh);
  mp=toastMPos(hMesh);
  patch('Vertices',vtx,'Faces',idx,'FaceVertexCData',col,'FaceColor', ...
        'interp','FaceAlpha',0.5,'EdgeColor',[0 0 0],'EdgeAlpha',0.5);
  axis equal;
  hold on
  for i=1:size(qp,1)
    vtx = sv;
    vtx(:,1) = vtx(:,1)+qp(i,1);
    vtx(:,2) = vtx(:,2)+qp(i,2);
    vtx(:,3) = vtx(:,3)+qp(i,3);
    patch('Vertices',vtx,'Faces',sf,'FaceVertexCData',sc,'FaceColor','interp','EdgeColor','none');
  end
  sc(:,2)=1;
  for i=1:size(mp,1)
    vtx = sv;
    vtx(:,1) = vtx(:,1)+mp(i,1);
    vtx(:,2) = vtx(:,2)+mp(i,2);
    vtx(:,3) = vtx(:,3)+mp(i,3);
    patch('Vertices',vtx,'Faces',sf,'FaceVertexCData',sc,'FaceColor','interp','EdgeColor','none');
  end
  %plot3(qp(:,1),qp(:,2),qp(:,3),'*r','MarkerSize',7)
  %plot3(mp(:,1),mp(:,2),mp(:,3),'ob','MarkerSize',7)

  daspect([1 1 1])
  view(3); axis tight
  camlight
  lighting gouraud
