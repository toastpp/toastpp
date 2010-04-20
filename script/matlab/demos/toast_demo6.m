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

% Last Modified by GUIDE v2.5 05-Jun-2008 12:24:39

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
function init(handles)

clear prm;
%prm.basis.meshfile = 'vox32_2blobs.msh';
%prm.basis.bdim = [25 30 20];
%prm.data.lnampfile = 'fmod_head_vox32_3plane.fem';
%prm.data.phasefile = 'farg_head_vox32_3plane.fem';
%prm.meas.qmfile = 'vox32_3plane.qm';
prm.fwdsolver.meshfile = 'cyl2.msh';
prm.solver.basis.bdim = [24 24 24];
prm.data.lnampfile = 'fmod_cyl2_5ring_100MHz.fem';
prm.data.phasefile = 'farg_cyl2_5ring_100MHz.fem';
prm.data.freq = 100;
prm.meas.qmfile = 'cyl_5ring.qm';
prm.meas.src = struct('type','Neumann','prof','Gaussian','width',2);
prm.meas.det = struct('prof','Gaussian','width',2);
prm.fwdsolver.method = 'BiCGSTAB';
%prm.linsolver.method = 'DIRECT';
prm.fwdsolver.tol = 1e-10;
prm.solver.method = 'PCG';
prm.solver.tol = 1e-8;
prm.solver.step0 = 87;
prm.solver.lsearch = true;
prm.regul.method = 'None';
prm.initprm.mua = struct('reset','HOMOG','val',0.01);
prm.initprm.mus = struct('reset','HOMOG','val',1);
prm.initprm.ref = struct('reset','HOMOG','val',1.4);

prm.callback.context = handles;
prm.callback.iter = @callback_vis; % iteration callback from recon

prm.basis.hMesh = toastReadMesh(prm.fwdsolver.meshfile);
toastReadQM (prm.basis.hMesh,prm.meas.qmfile);

prm.basis.hBasis = toastSetBasis('LINEAR',prm.basis.hMesh,prm.solver.basis.bdim,prm.solver.basis.bdim);
prm.smask = toastSolutionMask(prm.basis.hBasis);
n = toastMeshNodeCount(prm.basis.hMesh);
%load toast_demo6.mat
%prm.bmua = bmua; clear bmua;
%prm.bmus = bmus; clear bmus;
%prm.mua = toastMapBasisToMesh(prm.basis.hBasis,prm.bmua);
%prm.mus = toastMapBasisToMesh(prm.basis.hBasis,prm.bmus);
prm.mua = toastReadNIM('mua_tgt_cyl2.nim');
prm.mus = toastReadNIM('mus_tgt_cyl2.nim');
prm.ref = ones(n,1)*1.4;
prm.bmua = toastMapMeshToBasis(prm.basis.hBasis,prm.mua);
prm.bmus = toastMapMeshToBasis(prm.basis.hBasis,prm.mus);

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

showiso(handles,prm.bmua,0.015,prm.bmus,1.5,1);
showcs(handles,prm.bmua,lprm.cs_mua,0.005,0.02,handles.axes2);
showcs(handles,prm.bmus,lprm.cs_mus,0.5,2,handles.axes8);
tic

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if getappdata(handles.figure1,'isBusy') == false
    prm = getappdata(handles.figure1,'prm');
    setappdata(handles.figure1,'isBusy',true);
    set(handles.pushbutton1,'String','Stop (may be delayed)');
    drawnow
    toastRecon(prm);
    setappdata(handles.figure1,'isBusy',false);
    set(handles.pushbutton1,'String','Run');
else
    global ITR_TERM
    ITR_TERM = true;
end


% Display reconstruction results for current iteration
function callback_vis(handles,res)
prm = getappdata(handles.figure1,'prm');
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
semilogy(res.of); axis tight
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
switch get(hObject,'Value')
    case 1
        prm.regul = struct ('method','None');
    case 2
        prm.regul = struct ('method','TK1', ...
            'tau',1e-4, ...
            'prior',struct ('refname','','smooth',1,'threshold',0.1));
    case 3
        prm.regul = struct ('method','TV', ...
            'tau',1e-4, ...
            'prior',struct ('refname','','smooth',1,'threshold',0.1), ...
            'tv',struct ('beta',0.01));
    case 4
        prm.regul = struct ('method','Huber', ...
            'tau',1e-4, ...
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
%set(gcf,'Renderer','OpenGL');
cla;

[vtx,idx,perm] = toastSurfData(prm.basis.hMesh);
bb = toastMeshBB(prm.basis.hMesh);
pmax = bb(1,:);
pmin = bb(2,:);
nvtx = size(vtx,1);
msize=pmax-pmin;
col = zeros(nvtx,3); col(:,2)=1;
patch('Vertices',vtx,'Faces',idx,'FaceVertexCData',col,'FaceColor', ...
        'interp','EdgeColor','none','FaceAlpha',0.3);
bmua = reshape(bmua,nx,ny,nz);
[f,v] = isosurface(bmua,thres_mua);
if length(v) > 0
    v = [v(:,2) v(:,1) v(:,3)];
    for i=1:size(v,1)
        v(i,1) = v(i,1)/nx*msize(1) + pmin(1);
        v(i,2) = v(i,2)/ny*msize(2) + pmin(2);
        v(i,3) = v(i,3)/nz*msize(3) + pmin(3);
    end
    p = patch('Faces',f,'Vertices',v);
    isonormals(bmua,p);
    set(p,'FaceColor','red','EdgeColor','none');
end

bmus = reshape(bmus,nx,ny,nz);
[f,v] = isosurface(bmus,thres_mus);
if length(v) > 0
    v = [v(:,2) v(:,1) v(:,3)];
    for i=1:size(v,1)
        v(i,1) = v(i,1)/nx*msize(1) + pmin(1);
        v(i,2) = v(i,2)/ny*msize(2) + pmin(2);
        v(i,3) = v(i,3)/nz*msize(3) + pmin(3);
    end
    p = patch('Faces',f,'Vertices',v);
    isonormals(bmus,p);
    set(p,'FaceColor','blue','EdgeColor','none');
end

daspect([1 1 1])
view(-57,40);
axis equal;axis tight
camlight
lighting gouraud

% ============================================================
function showcs (handles,bimg,plane,imin,imax,hax)

prm = getappdata(handles.figure1,'prm');
nx = prm.solver.basis.bdim(1); ny = prm.solver.basis.bdim(2); nz = prm.solver.basis.bdim(3);
axes(hax);

img = zeros(size(bimg));
img(prm.smask) = bimg(prm.smask);

img = reshape(img,nx,ny,nz);
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
showcs(handles,prm.bmua,lprm.cs_mua,0.005,0.02,handles.axes2);
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
showcs(handles,prm.bmus,lprm.cs_mus,0.5,2,handles.axes8);
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


