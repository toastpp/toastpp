% N.B. To run this code you need to install Stanford Wavelab and add it to your matlab path
% (http://www-stat.stanford.edu/~wavelab/)
%
function FDOT_recon_CW_Proj

close all

nproj=8;
projgrid = [128 128];
meshfile = 'cyl2.msh';
qmfile = 'circle200_8x8x3_z0.qm';
mesh = toastMesh(meshfile);
mesh.ReadQM(qmfile);
hprojlist = toastMakeProjectorList(mesh,projgrid,[0 0 200],'ORTHO','pixelsize',0.25,'shift',[0 0]);
nh = mesh.NodeCount;
hpq = hprojlist(1:nproj);
grd = [32 32 32]
basis = toastBasis(mesh,grd);


%% Create the source vectors
qvec = mesh.Qvec ('Neumann', 'Gaussian', 2);


%% perform homogeneous forward calculations - Create the FEM system matrix

mua = ones(nh,1)*0.01;
mus = ones(nh,1)*1;
ref = ones(nh,1)*1.4;
smat = dotSysmat (mesh,mua,mus,ref,100);
smat = real(smat);

phi = full(smat\qvec);

figure;

for i=1:nproj
mesh.Display(log(squeeze(phi(:,i))));
pause(0.1)
end



%% fluorescence forward calculations

[vtx idx] = mesh.Data;
bb = mesh.BoundingBox;
bbmin = bb(1,:);
bbmax = bb(2,:);
cnt = [bbmax(1)/2-1 0 0]; rad = bbmax(2)/10;
fluo=zeros(nh,1);

for i=1:nh
    vtx(i,3)=0;
    if norm(vtx(i,1:3)-cnt)<=rad
        fluo(i) = 2;
    end
end



mfluo =repmat(fluo,1,nproj);

% for i=1:nproj
% qfluo(:,i) = toastIntFG(hmesh,mfluo(:,i)',phi(:,i)');
% end

qfluo = mfluo.*phi;
fl = smat\qfluo;

figure;

for i=1:nproj
mesh.Display(fl(:,i));
pause(0.1)
end



%% generate data
proj_ex = zeros([nproj projgrid]);

for k=1:nproj
    ex(:,k) = toastProjectToImage(hpq(k),phi(:,k));   % FLUORESCENCE DATA CCD Camera
    proj_ex(k,:,:) = reshape(ex(:,k),projgrid);
end

figure;
for i=1:nproj
imagesc(squeeze(proj_ex(i,:,:)));
pause(0.1)
end

proj_fluo = zeros([nproj projgrid]);

for k=1:nproj
    fl2(:,k) = toastProjectToImage(hpq(k),fl(:,k));   % FLUORESCENCE DATA CCD Camera
    proj_fluo(k,:,:) =reshape(fl2(:,k),projgrid);
end

figure;
for i=1:nproj
imagesc(squeeze(proj_fluo(i,:,:)));
pause(0.1)
end

%% Born Normalisation
 
immask = zeros([nproj projgrid]);
M = zeros([nproj projgrid]);
for k = 1:nproj
    exim = proj_ex(k,:,:); exim = exim(:);
    dmask = find( (exim > 0.1*max(exim)));
    im2 = zeros(projgrid(1:2));
    im2(dmask) = 1;
    M(k,:,:) = im2;
    im2(dmask) = im2(dmask)./exim(dmask);
    immask(k,:,:) = im2;
    subplot(1,3,1);imagesc(squeeze(proj_fluo(k,:,:)));colormap('jet');colorbar('horiz');title('data'); axis image;
    subplot(1,3,2);imagesc(squeeze(proj_fluo(k,:,:).*immask(k,:,:)));colorbar('horiz');title('mask'); axis image;
    subplot(1,3,3); imagesc(squeeze(proj_ex(k,:,:))); colorbar('horiz');title('excitation'); axis image;
    pause(1);
end

%% Data compression
dwidx = zeros(projgrid(1)*projgrid(2), nproj);
datawav = zeros(projgrid(1), projgrid(2), nproj);


qmf = MakeONFilter('Battle',3); %Wavelet filter
nw=10;


for q=1:nproj;
    datawav(:,:,q) = FWT2_PO(squeeze(proj_fluo(q,:,:).*immask(q,:,:)), 3, qmf);
    mask(:,:,q) = FWT2_PO(squeeze(immask(q,:,:)), 3, qmf);
    % sort data wavelet coeffs
    [dwsort, dwidx(:,q)] = sort(reshape(abs(datawav(:,:,q)),[],1));
    
end;

%--------------------------------------------------------------------------
%% JACOBIAN
%--------------------------------------------------------------------------

%% adjoint fields
data = zeros(nw*nproj,1);
J = zeros(nw*nproj, nh);
idx = 1;

for q=1:nproj
    display(['Calculating Jacobian: source ', num2str(q)])
    for k = size(dwidx,1)-nw+1:size(dwidx,1);
        mwav = zeros(projgrid);
        mwav(dwidx(k,q)) = 1;
        mwav_i = IWT2_PO(mwav, 3 , qmf);
        mvec2 = zeros([projgrid nproj]);
        mvec2(:,:,q) = mwav_i.*squeeze(immask(q,:,:));
        mv=toastProjectToField(hpq(q),mvec2(:,:,q));
        aphi = full(smat\mv);
        J(idx,:) = (phi(:,q)'.* aphi');
        %A=toastIntFG(hmesh,fluo',aphi');
        A=fluo'.*aphi';
        dy2(idx,:) =reshape((phi(:,q)'*(A')),[],1);
        dy(idx,:) = mv.'*fl(:,q);
        datawav_q = datawav(:,:,q);
        data(idx) = datawav_q(dwidx(k,q));   
        idx=idx+1;
    end
end

mfluo2= repmat(fluo,1,idx-1);
plot(diag((J'.'*mfluo2'')))
%% Reconstruction

H=J*J';
alpha =0.01*trace(H);
dlen = length(data);
x=J'*((H +alpha*speye(dlen))\(data(:)));
f = basis.Map('M->B',x);
%x=J'*((H +alpha*speye(dlen))\(dy(:)/cm));
%f = basis.Map('M->B',x);
fg = reshape(f,grd);
figure;
imagesc(fg(:,:,16));axis equal tight; colorbar
title('Fluo recon');

n = round(sqrt(grd(3)));

for i=1:grd(3)
subplot(n,n,i);
imagesc(fg(:,:,i));axis equal tight; 
end

