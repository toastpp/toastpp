function varargout = toast_demo2(varargin)
% TOAST_DEMO2 M-file for toast_demo2.fig
%      TOAST_DEMO2, by itself, creates a new TOAST_DEMO2 or raises the existing
%      singleton*.
%
%      H = TOAST_DEMO2 returns the handle to a new TOAST_DEMO2 or the handle to
%      the existing singleton*.
%
%      TOAST_DEMO2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TOAST_DEMO2.M with the given input arguments.
%
%      TOAST_DEMO2('Property','Value',...) creates a new TOAST_DEMO2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before toast_demo2_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to toast_demo2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help toast_demo2

% Last Modified by GUIDE v2.5 23-Feb-2008 04:19:30

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @toast_demo2_OpeningFcn, ...
                   'gui_OutputFcn',  @toast_demo2_OutputFcn, ...
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


% --- Executes just before toast_demo2 is made visible.
function toast_demo2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to toast_demo2 (see VARARGIN)

% Choose default command line output for toast_demo2
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

init(handles);
fwdsolve(handles);

% UIWAIT makes toast_demo2 wait for user response (see UIRESUME)
%uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = toast_demo2_OutputFcn(hObject, eventdata, handles) 
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
prm.basis.hMesh.ReadQM ('circle25_1x1.qm');

% initial source and detector positions
prm.phiq = 0;
prm.phim = pi;
radq = 24.5;
radm = 25.0;
Q = [radq*cos(prm.phiq) radq*sin(prm.phiq)];
M = [radm*cos(prm.phim) radm*sin(prm.phim)];
prm.basis.hMesh.SetQM(Q,M);

prm.qvec = prm.basis.hMesh.Qvec('Neumann','Gaussian',2);
prm.mvec = prm.basis.hMesh.Mvec('Gaussian',2,1.4);
prm.basis.hBasis = toastBasis(prm.basis.hMesh,[prm.bx prm.by],'Linear');
n = prm.basis.hMesh.NodeCount();
prm.mua = ones(n,1)*0.025;
prm.mus = ones(n,1)*2;
prm.ref = ones(n,1)*1.4;
prm.bkg = 0;
prm.meas.src.freq = 100;
prm.bmua = prm.basis.hBasis.Map('M->B',prm.mua);
prm.bmus = prm.basis.hBasis.Map('M->B',prm.mus);
axes(handles.axes1);
%h = surface(reshape(prm.bmua,prm.bx,prm.by),'EdgeColor','none'); axis tight
h = imagesc(rot90(reshape(prm.bmua,prm.bx,prm.by))); axis xy equal tight off
set(h,'ButtonDownFcn','toast_demo2(''axes1_ButtonDownFcn'',gcbo,[],guidata(gcbo))');
axes(handles.axes2);
h = imagesc(rot90(reshape(prm.bmus,prm.bx,prm.by))); axis xy equal tight off
%h = surface(reshape(prm.bmus,prm.bx,prm.by),'EdgeColor','none'); axis tight
set(h,'ButtonDownFcn','toast_demo2(''axes2_ButtonDownFcn'',gcbo,[],guidata(gcbo))');
axes(handles.axes7);
tmp = prm.basis.hBasis.Map('M->B',ones(n,1));
h = imagesc(rot90(reshape(tmp,prm.bx,prm.by))); axis xy equal tight off
%h = surface(reshape(tmp,prm.bx,prm.by),'EdgeColor','none'); axis tight
hold on;
x = prm.bx*(0.5+cos(prm.phiq)*0.45);
y = prm.by*(0.5-sin(prm.phiq)*0.45);
plot(x, y, 'og', 'MarkerSize',7, 'MarkerFaceColor','green');
x = prm.bx*(0.5+cos(prm.phim)*0.45);
y = prm.by*(0.5-sin(prm.phim)*0.45);
plot(x, y, 'oy', 'MarkerSize',7, 'MarkerFaceColor','yellow');
set(h,'ButtonDownFcn','toast_demo2(''axes7_ButtonDownFcn'',gcbo,[],guidata(gcbo))');
clear tmp

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
            load toast_demo2.mat
            prm.bmua = bmua; clear bmua;
            prm.bmus = bmus; clear bmus;
            prm.mua = prm.basis.hBasis.Map('B->M',prm.bmua);
            prm.mus = prm.basis.hBasis.Map('B->M',prm.bmus);
    end
    
    setappdata(handles.figure1,'prm',prm);
    fwdsolve(handles);
end


%% ===========================================================
function set_qm(handles,phiq,phim)
prm = getappdata(handles.figure1,'prm');

% Redefine QM positions
radq = 24.5; radm = 25.0;
Q = [radq*cos(phiq) radq*sin(phiq)];
M = [radm*cos(phim) radm*sin(phim)];
prm.basis.hMesh.SetQM(Q,M);
prm.qvec = prm.basis.hMesh.Qvec('Neumann','Gaussian',2);
prm.mvec = prm.basis.hMesh.Mvec('Gaussian',2,1.4);
prm.phiq = phiq;
prm.phim = phim;

axes(handles.axes7);
n = prm.basis.hMesh.NodeCount();
tmp = prm.basis.hBasis.Map('M->B',ones(n,1));
cla;h = imagesc(rot90(reshape(tmp,prm.bx,prm.by))); axis xy equal tight off
%cla;h = surface(reshape(tmp,prm.bx,prm.by),'EdgeColor','none'); axis tight
x = prm.bx*(0.5+cos(phiq)*0.45);
y = prm.by*(0.5-sin(phiq)*0.45);
hold on;
plot(x, y, 'og', 'MarkerSize',7, 'MarkerFaceColor','green');
x = prm.bx*(0.5+cos(phim)*0.45);
y = prm.by*(0.5-sin(phim)*0.45);
plot(x, y, 'oy', 'MarkerSize',7, 'MarkerFaceColor','yellow');
set(h,'ButtonDownFcn','toast_demo2(''axes7_ButtonDownFcn'',gcbo,[],guidata(gcbo))');
clear tmp

setappdata(handles.figure1,'prm',prm);
fwdsolve(handles);


%% ===========================================================
function setfreq(handles,freq)
prm = getappdata(handles.figure1,'prm');
if freq ~= prm.meas.src.freq
    prm.meas.src.freq = freq;
    set(handles.text20,'String',num2str(freq,'%0.0f'));

    setappdata(handles.figure1,'prm',prm);
    fwdsolve(handles);
end


%% ===========================================================
function fwdsolve(handles)
prm = getappdata(handles.figure1,'prm');

% perturbation parameters
pert = getappdata(handles.figure1,'pert');
bpert = zeros(prm.bx,prm.by);
bpert(pert.mua.x0-pert.mua.dx:pert.mua.x0+pert.mua.dx,pert.mua.y0-pert.mua.dx:pert.mua.y0+pert.mua.dx) = pert.mua.val;
bpert = reshape(bpert,[],1);
bmua = prm.bmua + bpert;
mua = prm.basis.hBasis.Map('B->M',bmua);
axes(handles.axes1);
cla;
h = imagesc(rot90(reshape(bmua,prm.bx,prm.by))); axis xy equal tight off
%h = surface(reshape(bmua,prm.bx,prm.by)','EdgeColor','none');axis tight
set(h,'ButtonDownFcn','toast_demo2(''axes1_ButtonDownFcn'',gcbo,[],guidata(gcbo))');

bpert = zeros(prm.bx,prm.by);
bpert(pert.mus.x0-pert.mus.dx:pert.mus.x0+pert.mus.dx,pert.mus.y0-pert.mus.dx:pert.mus.y0+pert.mus.dx) = pert.mus.val;
bpert = reshape(bpert,[],1);
bmus = prm.bmus + bpert;
mus = prm.basis.hBasis.Map('B->M',bmus);
axes(handles.axes2);
cla;
h = imagesc(rot90(reshape(bmus,prm.bx,prm.by))); axis xy equal tight off
%h = surface(reshape(bmus,prm.bx,prm.by)','EdgeColor','none');axis tight
set(h,'ButtonDownFcn','toast_demo2(''axes2_ButtonDownFcn'',gcbo,[],guidata(gcbo))');

% calculate the PMDFs
J = toastJacobian(prm.basis.hMesh, prm.basis.hBasis, prm.qvec, prm.mvec, ...
    mua, mus, prm.ref, prm.meas.src.freq, 'direct');
slen = size(J,2)/2;
Ja_lnamp = prm.basis.hBasis.Map('S->B', J(1,1:slen));
Jk_lnamp = prm.basis.hBasis.Map('S->B', J(1,slen+1:slen*2));
Ja_phase = prm.basis.hBasis.Map('S->B', J(2,1:slen));
Jk_phase = prm.basis.hBasis.Map('S->B', J(2,slen+1:slen*2));
axes(handles.axes3);
cla;imagesc(rot90(reshape(Ja_lnamp,prm.bx,prm.by))); axis xy equal tight off
%cla;surface(reshape(Ja_lnamp,prm.bx,prm.by)','EdgeColor','none'); axis tight;
axes(handles.axes4);
cla;imagesc(rot90(reshape(Jk_lnamp,prm.bx,prm.by))); axis xy equal tight off
%cla;surface(reshape(Jk_lnamp,prm.bx,prm.by)','EdgeColor','none'); axis tight;
axes(handles.axes5);
cla;imagesc(rot90(reshape(Ja_phase,prm.bx,prm.by))); axis xy equal tight off
%cla;surface(reshape(Ja_phase,prm.bx,prm.by)','EdgeColor','none'); axis tight;
axes(handles.axes6);
cla;imagesc(rot90(reshape(Jk_phase,prm.bx,prm.by))); axis xy equal tight off
%cla;surface(reshape(Jk_phase,prm.bx,prm.by)','EdgeColor','none'); axis tight;

function axes1_ButtonDownFcn(hObject, eventdata, handles)
prm = getappdata(handles.figure1,'prm');
pert = getappdata(handles.figure1,'pert');
mouse = get(gca,'currentpoint');
x = round(mouse(1,1));
y = prm.by-round(mouse(1,2));
x = min(max(x,pert.mua.dx+1),prm.bx-pert.mua.dx-1);
y = min(max(y,pert.mua.dx+1),prm.by-pert.mua.dx-1);
if x ~= pert.mua.x0 || y ~= pert.mua.y0
    pert.mua.x0 = x;
    pert.mua.y0 = y;
    setappdata(handles.figure1,'pert',pert);
    fwdsolve(handles);
end


function axes2_ButtonDownFcn(hObject, eventdata, handles)
prm = getappdata(handles.figure1,'prm');
pert = getappdata(handles.figure1,'pert');
mouse = get(gca,'currentpoint');
x = round(mouse(1,1));
y = prm.by-round(mouse(1,2));
x = min(max(x,pert.mus.dx+1),prm.bx-pert.mus.dx-1);
y = min(max(y,pert.mus.dx+1),prm.by-pert.mus.dx-1);
if x ~= pert.mus.x0 || y ~= pert.mus.y0
    pert.mus.x0 = x;
    pert.mus.y0 = y;
    setappdata(handles.figure1,'pert',pert);
    fwdsolve(handles);
end


function axes7_ButtonDownFcn(hObject, eventdata, handles)
prm = getappdata(handles.figure1,'prm');
mouse = get(gca,'currentpoint');
x = mouse(1,1);
y = prm.by-mouse(1,2);
dx = x/(prm.bx*0.5)-1.0;
dy = y/(prm.by*0.5)-1.0;
qx = cos(prm.phiq); qy = sin(prm.phiq);
mx = cos(prm.phim); my = sin(prm.phim);
dq = sqrt((qx-dx)^2 + (qy-dy)^2);
dm = sqrt((mx-dx)^2 + (my-dy)^2);
phi = atan2(dy,dx);
if dq < dm
    set_qm(handles,phi,prm.phim);
else
    set_qm(handles,prm.phiq,phi);
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
