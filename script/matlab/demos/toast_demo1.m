function varargout = toast_demo1(varargin)
% TOAST_DEMO1 M-file for toast_demo1.fig
%      TOAST_DEMO1, by itself, creates a new TOAST_DEMO1 or raises the
%      existing
%      singleton*.
%
%      H = TOAST_DEMO1 returns the handle to a new TOAST_DEMO1 or the handle to
%      the existing singleton*.
%
%      TOAST_DEMO1('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TOAST_DEMO1.M with the given input arguments.
%
%      TOAST_DEMO1('Property','Value',...) creates a new TOAST_DEMO1 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before toast_demo1_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to toast_demo1_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help toast_demo1

% Last Modified by GUIDE v2.5 23-Feb-2008 04:19:30

global FWD_IN_PROGRESS SRC_IN_PROGRESS;
FWD_IN_PROGRESS = false;
SRC_IN_PROGRESS = false;

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @toast_demo1_OpeningFcn, ...
                   'gui_OutputFcn',  @toast_demo1_OutputFcn, ...
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


% --- Executes just before toast_demo1 is made visible.
function toast_demo1_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to toast_demo1 (see VARARGIN)

% Choose default command line output for toast_demo1
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

init(handles);
fwdsolve(handles);

% UIWAIT makes toast_demo1 wait for user response (see UIRESUME)
%uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = toast_demo1_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


%% ===========================================================
function init(handles)
prm.bx = 64; prm.by = 64;
prm.blen = prm.bx*prm.by;
prm.basis.hMesh = toastMesh('circle25_32.msh');
prm.basis.hMesh.ReadQM('circle25_1x32.qm');
prm.basis.hBasis = toastBasis(prm.basis.hMesh,[prm.bx prm.by],'Linear');
prm.qvec = prm.basis.hMesh.Qvec('Neumann','Gaussian',2);
prm.mvec = prm.basis.hMesh.Mvec('Gaussian',2,1.4);
n = prm.basis.hMesh.NodeCount();
prm.mua = ones(n,1)*0.025;
prm.mus = ones(n,1)*2;
prm.ref = ones(n,1)*1.4;
prm.bkg = 0;
prm.meas.src.freq = 100;
prm.bmua = prm.basis.hBasis.Map('M->B',prm.mua);
prm.bmus = prm.basis.hBasis.Map('M->B',prm.mus);
prm.diff_fields = false;
prm.diff_proj = false;
axes(handles.axes1);
h = imagesc(rot90(reshape(prm.bmua,prm.bx,prm.by))); axis xy equal tight off
set(h,'ButtonDownFcn','toast_demo1(''axes1_ButtonDownFcn'',gcbo,[],guidata(gcbo))');
axes(handles.axes2);
h = imagesc(rot90(reshape(prm.bmus,prm.bx,prm.by))); axis xy equal tight off
set(h,'ButtonDownFcn','toast_demo1(''axes2_ButtonDownFcn'',gcbo,[],guidata(gcbo))');
axes(handles.axes7);
tmp = prm.basis.hBasis.Map('M->B',ones(n,1));
h = imagesc(rot90(reshape(tmp,prm.bx,prm.by))); axis xy equal tight off
x = prm.bx*(0.45+0.5);
y = prm.by*(0.5);
hold on; plot(x, y, 'og', 'MarkerSize',7, 'MarkerFaceColor','green');
set(h,'ButtonDownFcn','toast_demo1(''axes7_ButtonDownFcn'',gcbo,[],guidata(gcbo))');
clear tmp

% calculate unperturbed fields
smat = dotSysmat(prm.basis.hMesh,prm.mua,prm.mus,prm.ref,prm.meas.src.freq);
phi = smat\prm.qvec;
lphi = log(phi);
prm.lnamp0 = real(lphi);
prm.phase0 = -imag(lphi);

% calculate unperturbed projections
prm.lgamma0 = reshape (log(prm.mvec.' * phi), [], 1);

setappdata(handles.figure1,'prm',prm);

% perturbation parameters
pert.mua.x0 = round(prm.bx/4);
pert.mua.y0 = round(prm.by/2);
pert.mua.dx = 4;
pert.mua.val = 0.1;
pert.mus.x0 = round(prm.bx*3/4);
pert.mus.y0 = round(prm.by/2);
pert.mus.dx = 4;
pert.mus.val = 5;
setappdata(handles.figure1,'pert',pert);


%% ===========================================================
function setbackground(handles,bkg)
prm = getappdata(handles.figure1,'prm');
if bkg ~= prm.bkg
    prm.bkg = bkg;
    
    % update parameter background
    switch bkg
        case 0 % flat
            n = prm.basis.hMesh.NodeCount();
            prm.mua = ones(n,1)*0.025;
            prm.mus = ones(n,1)*2;
            prm.bmua = prm.basis.hBasis.Map('M->B',prm.mua);
            prm.bmus = prm.basis.hBasis.Map('M->B',prm.mus);
        case 1 % 'blobs'
            load toast_demo1.mat
            prm.bmua = bmua; clear bmua;
            prm.bmus = bmus; clear bmus;
            prm.mua = prm.basis.hBasis.Map('B->M',prm.bmua);
            prm.mus = prm.basis.hBasis.Map('B->M',prm.bmus);
    end
    
    % calculate unperturbed fields
    smat = dotSysmat(prm.basis.hMesh,prm.mua,prm.mus,prm.ref,...
                       prm.meas.src.freq);
    phi = smat\prm.qvec;
    lphi = log(phi);
    prm.lnamp0 = real(lphi);
    prm.phase0 = -imag(lphi);

    % calculate unperturbed projections
    prm.lgamma0 = reshape (log(prm.mvec.' * phi), [], 1);

    setappdata(handles.figure1,'prm',prm);
    fwdsolve(handles);
end


%% ===========================================================
function setsource(handles,phi)
global SRC_IN_PROGRESS;
SRC_IN_PROGRESS = true;
prm = getappdata(handles.figure1,'prm');

% Redefine Q position
M = prm.basis.hMesh.Mpos();
rad = 24.5;
Q = [rad*cos(phi) rad*sin(phi)];
prm.basis.hMesh.SetQM(Q,M);
prm.qvec = prm.basis.hMesh.Qvec('Neumann','Gaussian',2);
prm.mvec = prm.basis.hMesh.Mvec('Gaussian',2,1.4);

axes(handles.axes7);
n = prm.basis.hMesh.NodeCount();
tmp = prm.basis.hBasis.Map('M->B',ones(n,1));
cla;h = imagesc(rot90(reshape(tmp,prm.bx,prm.by))); axis xy equal tight off
x = prm.bx*(0.5+cos(phi)*0.45);
y = prm.by*(0.5-sin(phi)*0.45);
hold on; plot(x, y, 'og', 'MarkerSize',7, 'MarkerFaceColor','green');
set(h,'ButtonDownFcn','toast_demo1(''axes7_ButtonDownFcn'',gcbo,[],guidata(gcbo))');
clear tmp

% calculate unperturbed fields
smat = dotSysmat(prm.basis.hMesh,prm.mua,prm.mus,prm.ref,prm.meas.src.freq);
phi = smat\prm.qvec;
lphi = log(phi);
prm.lnamp0 = real(lphi);
prm.phase0 = -imag(lphi);

% calculate unperturbed projections
prm.lgamma0 = reshape (log(prm.mvec.' * phi), [], 1);

setappdata(handles.figure1,'prm',prm);
fwdsolve(handles);
SRC_IN_PROGRESS = false;


%% ===========================================================
function setfreq(handles,freq)
prm = getappdata(handles.figure1,'prm');
if freq ~= prm.meas.src.freq
    prm.meas.src.freq = freq;
    set(handles.text20,'String',num2str(freq,'%0.0f'));

    % calculate unperturbed fields
    smat = dotSysmat(prm.basis.hMesh,prm.mua,prm.mus,prm.ref,...
                       prm.meas.src.freq);
    phi = smat\prm.qvec;
    lphi = log(phi);
    prm.lnamp0 = real(lphi);
    prm.phase0 = -imag(lphi);

    % calculate unperturbed projections
    prm.lgamma0 = reshape (log(prm.mvec.' * phi), [], 1);

    setappdata(handles.figure1,'prm',prm);
    fwdsolve(handles);
end


%% ===========================================================
function fwdsolve(handles)
prm = getappdata(handles.figure1,'prm');

global FWD_IN_PROGRESS
FWD_IN_PROGRESS = true;

% perturbation parameters
pert = getappdata(handles.figure1,'pert');
bpert = zeros(prm.bx,prm.by);
bpert(pert.mua.x0-pert.mua.dx:pert.mua.x0+pert.mua.dx,pert.mua.y0-pert.mua.dx:pert.mua.y0+pert.mua.dx) = pert.mua.val;
bpert = reshape(bpert,[],1);
bmua = prm.bmua + bpert;
mua = prm.basis.hBasis.Map('B->M',bmua);
axes(handles.axes1);
cla; h = imagesc(rot90(reshape(bmua,prm.bx,prm.by))); axis xy equal tight off
set(h,'ButtonDownFcn','toast_demo1(''axes1_ButtonDownFcn'',gcbo,[],guidata(gcbo))');

bpert = zeros(prm.bx,prm.by);
bpert(pert.mus.x0-pert.mus.dx:pert.mus.x0+pert.mus.dx,pert.mus.y0-pert.mus.dx:pert.mus.y0+pert.mus.dx) = pert.mus.val;
bpert = reshape(bpert,[],1);
bmus = prm.bmus + bpert;
mus = prm.basis.hBasis.Map('B->M',bmus);
axes(handles.axes2);
cla; h = imagesc(rot90(reshape(bmus,prm.bx,prm.by))); axis xy equal tight off
set(h,'ButtonDownFcn','toast_demo1(''axes2_ButtonDownFcn'',gcbo,[],guidata(gcbo))');

% calculate the fields
n = prm.basis.hMesh.NodeCount();
smat = dotSysmat(prm.basis.hMesh,mua,mus,prm.ref,prm.meas.src.freq);
phi = smat\prm.qvec;
lphi = log(phi);
lnamp = real(lphi);
phase = -imag(lphi);
if prm.diff_fields == true;
    lnamp = lnamp - prm.lnamp0;
    phase = phase - prm.phase0;
end
axes(handles.axes3);
bphi = prm.basis.hBasis.Map('M->B',lnamp);
cla;imagesc(rot90(reshape(bphi,prm.bx,prm.by))); axis xy equal tight off;
axes(handles.axes4);
bphi = prm.basis.hBasis.Map('M->B',phase);
cla;imagesc(rot90(reshape(bphi,prm.bx,prm.by))); axis xy equal tight off;

% calculate projections
lgamma = reshape (log(prm.mvec.' * phi), [], 1);
axes(handles.axes5);
if prm.diff_proj == true
    plot(real(lgamma)-real(prm.lgamma0),'r');axis tight
else
    plot(real(prm.lgamma0));hold on;plot(real(lgamma),'r');hold off;axis tight
end
axes(handles.axes6);
if prm.diff_proj == true
    plot(-imag(lgamma)+imag(prm.lgamma0),'r');axis tight
else
    plot(-imag(prm.lgamma0));hold on;plot(-imag(lgamma),'r');hold off;axis tight
end
FWD_IN_PROGRESS = false;


function axes1_ButtonDownFcn(hObject, eventdata, handles)
prm = getappdata(handles.figure1,'prm');
h = gcf;
props.by = prm.by;
props.handles = handles;
props.axes = gca;
setappdata(h,'TestGuiCallbacks',props);
set(h,'WindowButtonMotionFcn',{@mua_ButtonMotionFcn});
set(h,'WindowButtonUpFcn',{@ButtonUpFcn});


function mua_ButtonMotionFcn(hObject, eventdata)
global FWD_IN_PROGRESS
props = getappdata(hObject,'TestGuiCallbacks');
handles = props.handles;
prm = getappdata(handles.figure1,'prm');
pert = getappdata(handles.figure1,'pert');
mouse = get(props.axes,'currentpoint');
x = round(mouse(1,1));
y = props.by-round(mouse(1,2));
x = min(max(x,pert.mua.dx+1),prm.bx-pert.mua.dx-1);
y = min(max(y,pert.mua.dx+1),prm.by-pert.mua.dx-1);
if (x ~= pert.mua.x0 || y ~= pert.mua.y0) && FWD_IN_PROGRESS == false
    pert.mua.x0 = x;
    pert.mua.y0 = y;
    setappdata(handles.figure1,'pert',pert);
    fwdsolve(handles);
end


function ButtonUpFcn(hObject, eventdata)
set(hObject,'WindowButtonMotionFcn','');
set(hObject,'WindowButtonUpFcn','');


function axes2_ButtonDownFcn(hObject, eventdata, handles)
prm = getappdata(handles.figure1,'prm');
h = gcf;
props.by = prm.by;
props.handles = handles;
props.axes = gca;
setappdata(h,'TestGuiCallbacks',props);
set(h,'WindowButtonMotionFcn',{@mus_ButtonMotionFcn});
set(h,'WindowButtonUpFcn',{@ButtonUpFcn});


function mus_ButtonMotionFcn(hObject, eventdata)
global FWD_IN_PROGRESS
props = getappdata(hObject,'TestGuiCallbacks');
handles = props.handles;
prm = getappdata(handles.figure1,'prm');
pert = getappdata(handles.figure1,'pert');
mouse = get(props.axes,'currentpoint');
x = round(mouse(1,1));
y = props.by-round(mouse(1,2));
x = min(max(x,pert.mus.dx+1),prm.bx-pert.mus.dx-1);
y = min(max(y,pert.mus.dx+1),prm.by-pert.mus.dx-1);
if (x ~= pert.mus.x0 || y ~= pert.mus.y0) && FWD_IN_PROGRESS == false
    pert.mus.x0 = x;
    pert.mus.y0 = y;
    setappdata(handles.figure1,'pert',pert);
    fwdsolve(handles);
end


function axes7_ButtonDownFcn(hObject, eventdata, handles)
h = gcf;
props.handles = handles;
props.axes = gca;
setappdata(h,'TestGuiCallbacks',props);
set(h,'WindowButtonMotionFcn',{@src_ButtonMotionFcn});
set(h,'WindowButtonUpFcn',{@ButtonUpFcn});


function src_ButtonMotionFcn(hObject, eventdata)
global SRC_IN_PROGRESS
props = getappdata(hObject,'TestGuiCallbacks');
handles = props.handles;
prm = getappdata(handles.figure1,'prm');
mouse = get(props.axes,'currentpoint');
x = mouse(1,1);
y = prm.by-mouse(1,2);
dx = x/(prm.bx*0.5)-1.0;
dy = y/(prm.by*0.5)-1.0;
phi = atan2(dy,dx);
if SRC_IN_PROGRESS == false
    setsource(handles,phi);
end


% --- Executes on slider movement.
function slider3_Callback(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
pert = getappdata(handles.figure1,'pert');
prm.bx = 64; prm.by = 64; % cheat
dx = round(get(hObject,'Value'));
if dx ~= pert.mua.dx
    pert.mua.dx = dx;
    pert.mua.x0 = min(max(pert.mua.x0,pert.mua.dx+1),prm.bx-pert.mua.dx-1);
    pert.mua.y0 = min(max(pert.mua.y0,pert.mua.dx+1),prm.by-pert.mua.dx-1);
    setappdata(handles.figure1,'pert',pert);
    fwdsolve(handles);
end


% --- Executes during object creation, after setting all properties.
function slider3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
%pert = getappdata(handles.figure1,'pert');
set(hObject,'Value',4);


% --- Executes on slider movement.
function slider4_Callback(hObject, eventdata, handles)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
pert = getappdata(handles.figure1,'pert');
val = get(hObject,'Value');
if val ~= pert.mua.val
    pert.mua.val = val;
    setappdata(handles.figure1,'pert',pert);
    fwdsolve(handles);
end


% --- Executes during object creation, after setting all properties.
function slider4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
set(hObject,'Value',0.1);


% --- Executes on slider movement.
function slider5_Callback(hObject, eventdata, handles)
% hObject    handle to slider5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
pert = getappdata(handles.figure1,'pert');
prm.bx = 64; prm.by = 64; % cheat
dx = round(get(hObject,'Value'));
if dx ~= pert.mus.dx
    pert.mus.dx = dx;
    pert.mus.x0 = min(max(pert.mus.x0,pert.mus.dx+1),prm.bx-pert.mus.dx-1);
    pert.mus.y0 = min(max(pert.mus.y0,pert.mus.dx+1),prm.by-pert.mus.dx-1);
    setappdata(handles.figure1,'pert',pert);
    fwdsolve(handles);
end


% --- Executes during object creation, after setting all properties.
function slider5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
set(hObject,'Value',4);


% --- Executes on slider movement.
function slider6_Callback(hObject, eventdata, handles)
% hObject    handle to slider6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
pert = getappdata(handles.figure1,'pert');
val = get(hObject,'Value');
if val ~= pert.mus.val
    pert.mus.val = val;
    setappdata(handles.figure1,'pert',pert);
    fwdsolve(handles);
end


% --- Executes during object creation, after setting all properties.
function slider6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
set(hObject,'Value',5);



% --- Executes on button press in radiobutton1.
function radiobutton1_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton1
prm = getappdata(handles.figure1,'prm');
if prm.diff_fields == true
    prm.diff_fields = false;
    set(handles.radiobutton2,'Value',0);
    setappdata(handles.figure1,'prm',prm);
    fwdsolve(handles);
end


% --- Executes on button press in radiobutton2.
function radiobutton2_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton2
prm = getappdata(handles.figure1,'prm');
if prm.diff_fields == false
    prm.diff_fields = true;
    set(handles.radiobutton1,'Value',0);
    setappdata(handles.figure1,'prm',prm);
    fwdsolve(handles);
end


% --- Executes on button press in radiobutton3.
function radiobutton3_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton3
prm = getappdata(handles.figure1,'prm');
if prm.diff_proj == true
    prm.diff_proj = false;
    set(handles.radiobutton4,'Value',0);
    setappdata(handles.figure1,'prm',prm);
    fwdsolve(handles);
end


% --- Executes on button press in radiobutton4.
function radiobutton4_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton4
prm = getappdata(handles.figure1,'prm');
if prm.diff_proj == false
    prm.diff_proj = true;
    set(handles.radiobutton3,'Value',0);
    setappdata(handles.figure1,'prm',prm);
    fwdsolve(handles);
end


% --- Executes on button press in radiobutton5.
function radiobutton5_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton5
prm = getappdata(handles.figure1,'prm');
if prm.bkg ~= 0
    set(handles.radiobutton6,'Value',0);
    setbackground(handles,0);
end


% --- Executes on button press in radiobutton6.
function radiobutton6_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton6
prm = getappdata(handles.figure1,'prm');
if prm.bkg ~= 1
    set(handles.radiobutton5,'Value',0);
    setbackground(handles,1);
end


% --- Executes on slider movement.
function slider7_Callback(hObject, eventdata, handles)
% hObject    handle to slider7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
prm = getappdata(handles.figure1,'prm');
freq = get(hObject,'Value');
setfreq(handles,freq);


% --- Executes during object creation, after setting all properties.
function slider7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
set(hObject,'Value',100);
