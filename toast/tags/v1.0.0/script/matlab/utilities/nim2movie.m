function nim2movie(nimfile,meshfile,moviename,fps,nimg,dim,imgmin,imgmax)
%
hMesh = toastReadMesh(meshfile);
hBasis = toastSetBasis(hMesh,[dim dim], [dim dim]);
for i = 1:nimg
    disp(num2str(i));
    nim = toastReadNIM(nimfile,i);
    img = toastMapMeshToGrid(hBasis,nim);
    img = reshape(img,dim,dim)';
    imagesc(img,[imgmin imgmax]);
    axis equal; axis tight;
    colormap(hot);%gray);
    drawnow;
    fname(i) = getframe(gcf);
end

movie2avi(fname,moviename,'FPS',fps);
