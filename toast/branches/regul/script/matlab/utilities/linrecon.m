function linrecon(varargin)

% LINRECON linear reconstruction program
%   LINRECON asks for options using a menu. The selections are saved to an
%   ini file.
%
%   LINRECON(FILENAME) looks for the options in a ini file, in ascii text
%   format. This can be used in a script to automate a series of
%   reconstructions or following changes made to an existing ini file.
%
%   The program expects difference data, consisting of log(amplitude) and
%   phase. It is not yet clear what the best way to estimate the standard
%   deviation is - try different methods for different data.
%
%   You can select whether to calculate the Jacobian (and solution mask) or
%   read them from a file by selecting 0 or 1 in the jac field.
%   If you're running Windows, you'll need to read them from a file. Once
%   it works, it might be a lot faster to save the Jacobian as a mat file
%   and load that in rather than a text file.
%
%   You can decide whether to reconstruct using amplitude and phase,
%   amplitude alone or phase alone by selecting 1, 2 or 3 in the datatypes
%   field.
%
%   Currently, the program normalises the Jacobian by the mean optical
%   properties, performs a row normalisation, divides the Jacobian and the
%   data by the standard deviation, and corrects for the coupling
%   coefficents. This might need tweaking depending on your data.
%
%   Display should work for 2D and 3D data but be prepared to modify that
%   if you have particular data to display.
%
%   Adam, 2 Aug 2006



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preliminaries
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
disp(pwd);

if isempty(varargin) || exist(varargin{:},'file')~=2

   %
   % No ini file selected - use menu
   %

   disp('No ini file selected - select from menu');

   % Dialog box for information
   prompt  = {'Mesh file :',...
       'QM file :',...
       'Log amplitude data file:',...
       'Log amplitude reference file :',...
       'Phase data file :',...
       'Phase reference file :',...
       'Pixel basis:',...
       'Starting parameters [mua mus n lambda frequency]',...
       'Recalculate Jacobian [1/0]',...
       'Amp+phase [1], amp only [2], phase only [3]',...
       'Enter root filename for output files :' };
   menutitle   = 'LINEAR IMAGE RECONSTRUCTION';
   lines   = 1;

   def     = {'mesh.opt','mod.qm','moddata.amp','modref.amp','moddata.phase','modref.phase','16 24 20','0.01 1 1.4 0.01 100','1','1','test'};
   answer  = inputdlg(prompt,menutitle,lines,def)';
   fields  = {'mesh','qm','ampdata','ampref','phasedata','phaseref','basis','start','jac','datatypes','savename'};

   A = cell2struct(answer, fields, 2);
   clear prompt menutitle def answer fields lines
   saveinifile(A, [A.savename '.ini']);

else

   %
   % ini file selected - use it
   %

   inifilename=varargin{:};
   disp(['Using ' inifilename]);
   A=loadinifile(inifilename);
end

disp(A);


% process filenames etc
numpix=str2num(A.basis);

% determine homogeneous starting values for mua and mus.
seperators=find(diff(isspace(A.start))>0);
mua=str2num(strtrim(A.start(1:seperators(1))));
mus=str2num(strtrim(A.start(seperators(1)+1:seperators(2))));
refind=str2num(A.start(seperators(2)+1:seperators(3)));
lambda=str2num(A.start(seperators(3)+1:seperators(4)));
frequency=str2num(A.start(seperators(4)+1:end));

kappa=1/(3*(mua+mus));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Load and process data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load qm file
[sourcepos, measpos, linklist] = loadqmfile(A.qm);

% load and process data files. It's unclear what the best way to calculate
% standard deviation is - possibilities are:
% if std(I)=sqrt(I)     std(ln(I))=sqrt(ln(I))
%                                  ln(sqrt(I))
%                                  1/sqrt(I)
% remember that the reconstruction requires ln(I). Check that this is your
% datatype.

if str2num(A.datatypes)==1 || str2num(A.datatypes)==2
   absamp = loaddatafile(A.ampdata);
   refamp =  loaddatafile(A.ampref);
   amp = absamp-refamp;
   stdevamp=log(sqrt(exp(refamp)));
end

if str2num(A.datatypes)==1 || str2num(A.datatypes)==3
   absphase =  loaddatafile(A.phasedata);
   refphase =  loaddatafile(A.phaseref);
   phase = (absphase-refphase);
   stdevphase = abs(refphase*0.01);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate Jacobian using matlab toast or read from a file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if str2num(A.jac)
   hMesh = toastReadMesh (A.mesh);
   toastReadQM (hMesh,A.qm );

   % Set up the mapper between FEM and solution bases
   hBasis = toastSetBasis (hMesh, numpix, numpix*2);
   solmask = toastSolutionMask (hBasis);

   % Generate source and measurement vectors
   qvec = toastQvec (hMesh, 'Neumann', 'Gaussian', 2);
   mvec = toastMvec (hMesh, 'Gaussian', 2);

   % Construct the Jacobian
   tic
   dummy=ones(length(qvec),1);
   keyboard
   J = toastJacobian (hMesh, hBasis, qvec, mvec, mua*dummy, mus*dummy, refind*dummy, 100, 'gmres');
   clear mvec qvec dummy
   toc

   stats(J);
   save J.txt J -double -ascii -tabs
   save solmask.txt solmask -double -ascii -tabs
   toc
else
   load J.txt;
   load solmask.txt
end


switch A.datatypes
   case '1'
       disp('Reconstruct with amp and phase');
       data = [amp phase];
       stdev = [stdevamp stdevphase];
   case '2'
       disp('Reconstruct with amp only');
       data = amp;
       stdev = stdevamp;
       J=J(1:end/2,:);
   case '3'
       disp('Reconstruct with phase only');
       data = phase;
       stdev = stdevphase;
       J=J(end/2+1:end,:);
end

[nmeas, nnodes]=size(J);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% normalise Jacobian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% nomalise optical parameters - no effect if only one parameter
if str2num(A.datatypes)==1
   disp('normalising optical properties');
   J=J.*[repmat(ones(1, size(J,2)/2)*mua, nmeas, 1) repmat(ones(1, size(J,2)/2)*kappa, nmeas, 1)];
end


% row normalisation
disp('row normalisation');
R1=1./sqrt(sum(J.^2,2))';
J=[repmat(R1,nnodes,1)'].*J;

% stdev
disp('normalise by stdev');
R2=1./stdev;
J=[repmat(R2,nnodes,1)'].*J;
data=data.*R2;

% Append correction for coupling coefficients onto J
S=zeros(nmeas,size(sourcepos,1));
D=zeros(nmeas,size(measpos,1));
L=linklist+1;
measnum=0;
for ss=1:size(L,1)
   for dd=1:size(L,2)
       if  L(ss,dd)>0
           measnum=measnum+1;
           S(measnum,ss)=1;
           D(measnum,L(ss,dd))=1;
       end
   end
end

J=[J S D];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reconstruction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


disp('reconstruction...');
JJT=J*J';
figure;imagesc(JJT);
S=svd(JJT);
figure;semilogy(S/max(S));
lambda1=lambda*max(S);

invJ = J' * inv( JJT + eye(length(JJT)).*lambda1);
img = invJ * data';
imgupdate=img(1:nnodes,:);
CC=img([nnodes+1:end],:);


% convert image in mua and kappa to mua and mus
muaupdate = imgupdate(1:nnodes/2);
kappaupdate = imgupdate(nnodes/2+1:end);

muab = mua * exp (muaupdate/mua);
kappab = kappa * exp (kappaupdate/kappa);
musb = 1./(3*kappab)-mua;


disp(' ');
disp('image update');stats(imgupdate);
if length(CC)>0
   disp('coupling coefficients');stats(CC);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  DISPLAY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

muaimg=zeros(numpix);
musimg=zeros(numpix);
muaimg(solmask)=muab;
musimg(solmask)=musb;

map=colormap('jet');
%map=colormap('gray');

if ndims(muaimg)==2 % 2D images
   % mua image
   h1=figure;
   imgstat=muaimg(muaimg>0);
   maxmua=max(imgstat(:));
   minmua=min(imgstat(:));
   scale=0.5*length(map)/(max(abs(muaimg(:))-mua))
   slice = (scale*(muaimg-mua))+length(map)/2;
   image(slice');axis equal;axis tight
   ax = gca;
   axes('position',[.1 .86 .8 .05],'Box','off','Visible','off');
   title(['\mu_a: ' num2str(minmua) ' - ' num2str(maxmua) ' mm^{-1}'] );
   set(get(gca,'Title'),'Visible','On');
   set(get(gca,'Title'),'FontSize',18);
   axes(ax);
   saveas(gcf,[A.savename '_mua.png'],'png');

   % mus image
   h2=figure;
   imgstat=musimg(musimg>0);
   maxmus=max(imgstat(:));
   minmus=min(imgstat(:));
   scale=0.5*length(map)/(max(abs(musimg(:))-mus))
   slice = (scale*(musimg-mus))+length(map)/2;
   image(slice');axis equal;axis tight
   ax = gca;
   axes('position',[.1 .86 .8 .05],'Box','off','Visible','off');
   title(['\mu_a: ' num2str(minmus) ' - ' num2str(maxmus) ' mm^{-1}'] );
   set(get(gca,'Title'),'Visible','On');
   set(get(gca,'Title'),'FontSize',18);
   axes(ax);
   saveas(gcf,[A.savename '_mus.png'],'png');
end


if ndims(muaimg)==3 % 3d image

   dim=3; % change this for different views

   slices=1:2:size(muaimg,dim);
   nx=ceil(sqrt(length(slices)));
   ny=ceil(length(slices)/nx);

   % mua
   h1=figure;
   muaimg=shiftdim(muaimg,dim);
   IMG=muaimg(:,end:-1:1,slices);
   imgstat=IMG(IMG>0);
   maxmua=max(imgstat(:));
   minmua=min(imgstat(:));
   scale=0.5*length(map)/(max(abs(IMG(:))-mua))

   for sl=1:length(slices)
       slice = squeeze(IMG(:,:,sl));
       slice = (scale*(slice-mua))+length(map)/2;
       subplot(nx,ny,sl);image(slice');axis equal;axis tight
   end

   ax = gca;
   axes('position',[.1 .86 .8 .05],'Box','off','Visible','off');
   title(['\mu_a: ' num2str(minmua) ' - ' num2str(maxmua) ' mm^{-1}'] );
   set(get(gca,'Title'),'Visible','On');
   set(get(gca,'Title'),'FontSize',18);
   axes(ax);
   saveas(gcf,[A.savename '_mua.png'],'png');

   % mus
   h2=figure;
   musimg=shiftdim(musimg,dim);
   IMG=musimg(:,end:-1:1,slices);
   imgstat=IMG(IMG>0);
   maxmus=max(imgstat(:));
   minmus=min(imgstat(:));
   scale=0.5*length(map)/(max(abs(IMG(:))-mus))

   for sl=1:length(slices)
       slice = squeeze(IMG(:,:,sl));
       slice = (scale*(slice-mus))+length(map)/2;
       H(sl)=subplot(nx,ny,sl);image(slice');axis equal;axis tight
   end

   ax = gca;
   axes('position',[.1 .86 .8 .05],'Box','off','Visible','off');
   title(['\mu`_s: ' num2str(minmus) ' - ' num2str(maxmus) ' mm^{-1}'] );
   set(get(gca,'Title'),'Visible','On');
   set(get(gca,'Title'),'FontSize',18);
   axes(ax);
   saveas(gcf,[A.savename '_mus.png'],'png');
end

op=fopen([A.savename 'stats.txt'],'w');
fprintf(op,'mua %f - %f\n',minmua,maxmua);
fprintf(op,'mus %f - %f\n',minmus,maxmus);
fclose(op);
