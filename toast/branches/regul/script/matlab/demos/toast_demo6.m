function varargout = toast_demo6(varargin)
% TOAST_DEMO6 M-file for toast_demo6.fig
%      TOAST_DEMO6, by itself, creates a new TOAST_DEMO6 or raises the existing
%      singleton*.
%
%      H = TOAST_DEMO6 returns the handle to a new TOAST_DEMO6 or the handle to
%      the existing singleton*.
%
%      TOAST_DEMO6('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TOAST_DEMO6.M with the given input arguments.
%
%      TOAST_DEMO6('Property','Value',...) creates a new TOAST_DEMO6 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before toast_demo6_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to toast_demo6_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help toast_demo6

% Last Modified by GUIDE v2.5 27-Oct-2014 17:33:06

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @toast_demo6_OpeningFcn, ...
                   'gui_OutputFcn',  @toast_demo6_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before toast_demo6 is made visible.
function toast_demo6_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to toast_demo6 (see VARARGIN)

% Choose default command line output for toast_demo6
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

init(handles);

% UIWAIT makes toast_demo6 wait for user response (see UIRESUME)
%uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = toast_demo6_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


%% ===========================================================

function prm = load_target(prm)
meshdir = '../../../test/3D/meshes/';
meshfile = 'cyl4_blobs.msh';
qmfile = prm.meas.qmfile;
prm.user.target.mesh = toastMesh([meshdir meshfile]);
prm.user.target.mesh.ReadQM(qmfile);
prm.user.target.basis = toastBasis(prm.user.target.mesh, prm.solver.basis.bdim);
nim = toastNim([meshdir 'mua_tgt_cyl4.nim']);
prm.user.target.mua = nim.Values;
prm.user.target.bmua = prm.user.target.basis.Map('M->B', prm.user.target.mua);
nim = toastNim([meshdir 'mus_tgt_cyl4.nim']);
prm.user.target.mus = nim.Values;
prm.user.target.bmus = prm.user.target.basis.Map('M->B', prm.user.target.mus);

%% ===========================================================

function prm = gen_data(prm)
mesh = prm.user.target.mesh;
n = mesh.NodeCount;
mua = prm.user.target.mua;
mus = prm.user.target.mus;
ref = ones(n,1) * 1.4;
freq = prm.data.freq;
smat = dotSysmat(mesh,mua,mus,ref,freq);
qvec = mesh.Qvec(prm.meas.src.type,prm.meas.src.prof,prm.meas.src.width);
mvec = mesh.Mvec(prm.meas.det.prof,prm.meas.det.width,ref);
data = reshape(mvec.' * (smat\qvec), [], 1);
lndata = log(data);
% add noise
lndata = lndata + lndata.*prm.user.noiselevel.*randn(size(lndata));
prm.data.lnamp = real(lndata);
prm.data.phase = imag(lndata);

%% ===========================================================

function init(handles)

meshdir = '../../../test/3D/meshes/';
prm = toastParam;
%prm.basis.meshfile = 'vox32_2blobs.msh';
%prm.basis.bdim = [25 30 20];
%prm.data.lnampfile = 'fmod_head_vox32_3plane.fem';
%prm.data.phasefile = 'farg_head_vox32_3plane.fem';
%prm.meas.qmfile = 'vox32_3plane.qm';
prm.fwdsolver.meshfile = [meshdir 'cyl2.msh'];
prm.solver.basis.bdim = [32 32 32];
%prm.data.lnampfile = 'fmod_cyl2_5ring_100MHz.fem';
%prm.data.phasefile = 'farg_cyl2_5ring_100MHz.fem';
prm.data.freq = 100;
%prm.meas.qmfile = 'cyl_5ring.qm';
prm.meas.qmfile = [meshdir 'cyl_5ring.qm'];
prm.meas.src = struct('type','Neumann','prof','Gaussian','width',2);
prm.meas.det = struct('prof','Gaussian','width',2);
%prm.fwdsolver.method = 'BiCGSTAB';
prm.fwdsolver.method = 'DIRECT';
prm.fwdsolver.tol = 1e-10;
prm.fwdsolver.hmesh = toastMesh(prm.fwdsolver.meshfile);
prm.fwdsolver.hmesh.ReadQM (prm.meas.qmfile);
prm.solver.method = 'PCG';
prm.solver.tol = 1e-8;
prm.solver.step0 = 1e-2;
prm.solver.lsearch = true;
prm.regul.method = 'None';
prm.initprm.mua = struct('reset','HOMOG','val',0.01);
prm.initprm.mus = struct('reset','HOMOG','val',1);
prm.initprm.ref = struct('reset','HOMOG','val',1.4);

prm.transient.callback.context = handles;
prm.transient.callback.iter = @callback_vis; % iteration callback from recon

prm.solver.basis.hbasis = toastBasis(prm.fwdsolver.hmesh,prm.solver.basis.bdim,'Linear');
%prm.smask = toastSolutionMask(prm.basis.hBasis);
n = prm.fwdsolver.hmesh.NodeCount;

% load target data
prm = load_target(prm);

%prm.initprm.mua.mua = toastNim('mua_tgt_cyl2.nim');
%prm.initprm.mus.mus = toastNim('mus_tgt_cyl2.nim');
prm.initprm.mua.mua = ones(n,1)*0.01;
prm.initprm.mus.mus = ones(n,1)*1;
prm.initprm.ref.ref = ones(n,1)*1.4;
prm.initprm.mua.bmua = prm.solver.basis.hbasis.Map('M->B',prm.initprm.mua.mua);
prm.initprm.mus.bmus = prm.solver.basis.hbasis.Map('M->B',prm.initprm.mus.mus);

setappdata(handles.figure1,'prm',prm);

set(handles.pushbutton1,'String','Run');
setappdata(handles.figure1,'isBusy',false);

lprm.cs_mua = round(prm.solver.basis.bdim(3)*0.2333);
lprm.cs_mus = round(prm.solver.basis.bdim(3)*0.7333);
lprm.res_mua = [];
lprm.res_mus = [];
sstep = 1/(prm.solver.basis.bdim(3)-1);
set(handles.slider1,'Min',1);
set(handles.slider1,'Max',prm.solver.basis.bdim(3));
set(handles.slider1,'SliderStep',[sstep,sstep*5]);
set(handles.slider1,'Value',lprm.cs_mua);
set(handles.slider2,'Min',1);
set(handles.slider2,'Max',prm.solver.basis.bdim(3));
set(handles.slider2,'SliderStep',[sstep,sstep*5]);
set(handles.slider2,'Value',lprm.cs_mus);

setappdata(handles.figure1,'lprm',lprm);

showiso(handles,prm.user.target.bmua,0.015,prm.user.target.bmus,1.5,1);
showcs(handles,prm.user.target.bmua,lprm.cs_mua,0.005,0.02,handles.axes2);
showcs(handles,prm.user.target.bmus,lprm.cs_mus,0.5,2,handles.axes8);
tic

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if getappdata(handles.figure1,'isBusy') == false
    prm = getappdata(handles.figure1,'prm');
    prm.regul.tau = sscanf(get(handles.edit1,'String'),'%f');
    prm.user.noiselevel = sscanf(get(handles.edit2,'String'),'%f');
    setappdata(handles.figure1,'isBusy',true);
    set(handles.pushbutton1,'String','Stop (may be delayed)');
    drawnow
    % generate the forward data
    prm = gen_data(prm);
    toastRecon(prm);
    setappdata(handles.figure1,'isBusy',false);
    set(handles.pushbutton1,'String','Run');
else
    global ITR_TERM
    ITR_TERM = true;
end


% Display reconstruction results for current iteration
function callback_vis(prm,res)
handles = prm.transient.callback.context;
lprm = getappdata(handles.figure1,'lprm');
th_mua = min (0.015, max (0.011,(max(res.bmua)+0.01)/2));
th_mus = min (1.5, max (1.1, (max(res.bmus)+1)/2));
showiso(handles,res.bmua,th_mua,res.bmus,th_mus,3);
showcs(handles,res.bmua,lprm.cs_mua,max(0.005,min(res.bmua)),max(res.bmua),handles.axes4);
showcs(handles,res.bmus,lprm.cs_mus,max(0.5,min(res.bmus)),max(res.bmus),handles.axes9);
lprm.res_mua = res.bmua;
lprm.res_mus = res.bmus;
setappdata(handles.figure1,'lprm',lprm);

axes(handles.axes5);
plot(res.of); axis tight
set(handles.axes5,'FontSize',7);
toc
tic
drawnow


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1
prm = getappdata(handles.figure1,'prm');
method = {'PCG' 'LM'};
prm.solver.method = method{get(hObject,'Value')};
setappdata(handles.figure1,'prm',prm);


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2
prm = getappdata(handles.figure1,'prm');
switch get(handles.popupmenu2,'Value')
    case 1
        prm.regul = struct ('method','None');
    case 2
        prm.regul = struct ('method','TK1', ...
            'prior',struct ('refname','','smooth',1,'threshold',0.1));
    case 3
        prm.regul = struct ('method','TV', ...
            'prior',struct ('refname','','smooth',1,'threshold',0.1), ...
            'tv',struct ('beta',0.01));
    case 4
        prm.regul = struct ('method','Huber', ...
            'prior',struct ('refname','','smooth',1,'threshold',0.1), ...
            'huber',struct ('eps',0.01));
end
setappdata(handles.figure1,'prm',prm);


% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function showiso(handles,bmua,thres_mua,bmus,thres_mus,fignum)
  
prm = getappdata(handles.figure1,'prm');
nx = prm.solver.basis.bdim(1); ny = prm.solver.basis.bdim(2); nz = prm.solver.basis.bdim(3);
eval(['axes(handles.axes' num2str(fignum) ')']);
cla;
set(gcf,'Renderer','OpenGL');

[vtx,idx,perm] = prm.fwdsolver.hmesh.SurfData();
bb = prm.fwdsolver.hmesh.BoundingBox();
pmin = bb(1,:);
pmax = bb(2,:);
nvtx = size(vtx,1);
msize=pmax-pmin;
[x,y,z] = meshgrid([pmin(1):msize(1)/(nx-1):pmax(1)], ...
                   [pmin(2):msize(2)/(ny-1):pmax(2)], ...
                   [pmin(3):msize(3)/(nz-1):pmax(3)]);
               
bmua = reshape(bmua,nx,ny,nz);
[f,v] = isosurface(x,y,z,bmua,thres_mua);
if length(v) > 0
    p = patch('Faces',f,'Vertices',v);
    isonormals(x,y,z,bmua,p);
    set(p,'FaceColor',[1 0.5 0.5],'EdgeColor','none');
end

bmus = reshape(bmus,nx,ny,nz);
[f,v] = isosurface(x,y,z,bmus,thres_mus);
if length(v) > 0
    p = patch('Faces',f,'Vertices',v);
    isonormals(x,y,z,bmus,p);
    set(p,'FaceColor',[0.5 0.5 1],'EdgeColor','none');
end

col = zeros(nvtx,3); col(:,2)=1;
patch('Vertices',vtx,'Faces',idx,'FaceVertexCData',col,'FaceColor', ...
        'interp','EdgeColor','none','FaceAlpha',0.2);

daspect([1 1 1])
view(-57,40);
axis equal;axis tight
light('Position',[-20,0,0]);
lighting gouraud

% ============================================================
function showcs (handles,bimg,plane,imin,imax,hax)

prm = getappdata(handles.figure1,'prm');
nx = prm.solver.basis.bdim(1); ny = prm.solver.basis.bdim(2); nz = prm.solver.basis.bdim(3);
axes(hax);

%img = zeros(size(bimg));
%img(prm.smask) = bimg(prm.smask);

img = reshape(bimg,nx,ny,nz);
img = squeeze(img(:,:,plane));

imagesc(img',[imin,imax]);
axis('xy');
axis equal tight off
%set(hax,'XTick',[],'XTickLabel','','YTick',[],'YTickLabel','');


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
prm = getappdata(handles.figure1,'prm');
lprm = getappdata(handles.figure1,'lprm');
lprm.cs_mua = round(get(hObject,'Value'));
setappdata(handles.figure1,'lprm',lprm);
bmua = prm.user.target.bmua;
showcs(handles,bmua,lprm.cs_mua,0.005,0.02,handles.axes2);
if length(lprm.res_mua) > 0
    showcs(handles,lprm.res_mua,lprm.cs_mua,max(0.005,min(lprm.res_mua)),min(0.02,max(lprm.res_mua)),handles.axes4);
end


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider2_Callback(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
prm = getappdata(handles.figure1,'prm');
lprm = getappdata(handles.figure1,'lprm');
lprm.cs_mus = round(get(hObject,'Value'));
setappdata(handles.figure1,'lprm',lprm);
bmus = prm.user.target.bmus;
showcs(handles,bmus,lprm.cs_mus,0.5,2,handles.axes8);
if length(lprm.res_mus) > 0
    showcs(handles,lprm.res_mus,lprm.cs_mus,max(0.5,min(lprm.res_mus)),min(2,max(lprm.res_mus)),handles.axes9);
end


% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end





function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
