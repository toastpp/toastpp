function varargout = toast_demo5(varargin)
% TOAST_DEMO5 M-file for toast_demo5.fig
%      TOAST_DEMO5, by itself, creates a new TOAST_DEMO5 or raises the existing
%      singleton*.
%
%      H = TOAST_DEMO5 returns the handle to a new TOAST_DEMO5 or the handle to
%      the existing singleton*.
%
%      TOAST_DEMO5('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TOAST_DEMO5.M with the given input arguments.
%
%      TOAST_DEMO5('Property','Value',...) creates a new TOAST_DEMO5 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before toast_demo5_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to toast_demo5_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help toast_demo5

% Last Modified by GUIDE v2.5 07-Mar-2008 23:20:19

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @toast_demo5_OpeningFcn, ...
                   'gui_OutputFcn',  @toast_demo5_OutputFcn, ...
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


% --- Executes just before toast_demo5 is made visible.
function toast_demo5_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to toast_demo5 (see VARARGIN)

% Choose default command line output for toast_demo5
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

init(handles);

% UIWAIT makes toast_demo5 wait for user response (see UIRESUME)
%uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = toast_demo5_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


%% ===========================================================
function init(handles)

clear prm;
%prm.fwdsolver.meshfile = 'circle25_32.msh';
prm.fwdsolver.meshfile = 'ellips_tri3.msh';
prm.solver.basis.bdim = [96 96];
prm.data.lnampfile = 'fmod_ellips_16x16_100MHz.fem';
prm.data.phasefile = 'farg_ellips_16x16_100MHz.fem';
prm.data.freq = 100;
prm.meas.qmfile = 'circle25_16x16.qm';
prm.meas.src = struct('type','Neumann','prof','Gaussian','width',2);
prm.meas.det = struct('prof','Gaussian','width',2);
prm.fwdsolver.method = 'Direct';
prm.fwdsolver.tol = 1e-10;
prm.solver.method = 'LM';
prm.solver.tol = 1e-8;
prm.solver.step0 = 100;
load toast_demo5
prm.bmua = bmua_tgt;
prm.bmus = bmus_tgt;
prm.regul = struct ('method','TK1', ...
     'tau',1e-4, ...
     'prior',struct ('smooth',1,'threshold',0.1));
prm.regul.tv.beta = 0.01;
prm.regul.huber.eps = 0.01;
clear bmua_tgt bmus_tgt bkap_tgt;
prm.initprm.mua = struct('reset','HOMOG','val',0.025);
prm.initprm.mus = struct('reset','HOMOG','val',2);
prm.initprm.ref = struct('reset','HOMOG','val',1.4);

prm.callback.context = handles;
prm.callback.iter = @callback_vis; % iteration callback from recon
prm.callback.request.prior = true;

prm.basis.hMesh = toastReadMesh(prm.fwdsolver.meshfile);
toastReadQM (prm.basis.hMesh,prm.meas.qmfile);

prm.basis.hBasis = toastSetBasis('LINEAR',prm.basis.hMesh,prm.solver.basis.bdim,prm.solver.basis.bdim);
prm.smask = toastSolutionMask(prm.basis.hBasis);
n = toastMeshNodeCount(prm.basis.hMesh);
prm.mua = toastMapBasisToMesh(prm.basis.hBasis,prm.bmua);
prm.mus = toastMapBasisToMesh(prm.basis.hBasis,prm.bmus);
prm.ref = ones(n,1)*1.4;
axes(handles.axes1);
imagesc(reshape(prm.bmua,prm.solver.basis.bdim(1),prm.solver.basis.bdim(2)),[0.005 0.05]);
axis tight
set(handles.axes1,'XTick',[],'XTickLabel','','YTick',[],'YTickLabel','');
axes(handles.axes2);
imagesc(reshape(prm.bmus,prm.solver.basis.bdim(1),prm.solver.basis.bdim(2)),[0.5 4]);
axis tight
set(handles.axes2,'XTick',[],'XTickLabel','','YTick',[],'YTickLabel','');

setappdata(handles.figure1,'prm',prm);
showref(handles,prm);

set(handles.pushbutton1,'String','Run');
setappdata(handles.figure1,'isBusy',false);


%% --- Display the diffusivity reference images
function showref(handles,prm)
% create a regularisation instance on the fly to get the
% diffusivity prior
cm = 0.3/1.4;
slen = length(prm.smask);
scmua = ones(slen,1)*0.025*cm;
sckap = ones(slen,1)*(1/(3*2.025))*cm;
logx = log([scmua;sckap]);
hReg = toastRegul (prm.regul, prm.basis.hBasis, logx);
res.kapref = toastRegulKappa (hReg,logx);
toastDeleteRegul (hReg);
disp_kapref (handles,res,prm);


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if getappdata(handles.figure1,'isBusy') == false
    prm = getappdata(handles.figure1,'prm');
    setappdata(handles.figure1,'isBusy',true);
    set(handles.pushbutton1,'String','Stop (may be delayed)');
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

axes(handles.axes3);
imagesc(reshape(res.bmua,res.bdim(1),res.bdim(2)),[0.005 0.05]); axis tight
set(handles.axes3,'XTick',[],'XTickLabel','','YTick',[],'YTickLabel','');

axes(handles.axes4);
imagesc(reshape(res.bmus,res.bdim(1),res.bdim(2)),[0.5 4]); axis tight
set(handles.axes4,'XTick',[],'XTickLabel','','YTick',[],'YTickLabel','');

axes(handles.axes5);
semilogy(res.of); axis tight
set(handles.axes5,'FontSize',7);

disp_kapref(handles,res,prm);
drawnow


function disp_kapref(handles,res,prm)
slen = length(res.kapref)/2;

axes(handles.axes7);
kref = toastMapSolToBasis(prm.basis.hBasis,res.kapref(1:slen));
imagesc(reshape(kref,prm.solver.basis.bdim(1),prm.solver.basis.bdim(2)),[0 1.2]); axis tight
set(handles.axes7,'XTick',[],'XTickLabel','','YTick',[],'YTickLabel','');

axes(handles.axes8);
kref = toastMapSolToBasis(prm.basis.hBasis,res.kapref(slen+1:slen*2));
imagesc(reshape(kref,prm.solver.basis.bdim(1),prm.solver.basis.bdim(2)),[0 1.2]); axis tight
set(handles.axes8,'XTick',[],'XTickLabel','','YTick',[],'YTickLabel','');


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1

prm = getappdata(handles.figure1,'prm');

switch get(hObject,'Value')
    case 1
        if isfield(prm.regul.prior,'refimg')
            prm.regul.prior = rmfield(prm.regul.prior,'refimg');
        end
    case 2
        load toast_demo5 bmua_tgt bkap_tgt
        prm.regul.prior.refimg = [bmua_tgt; bkap_tgt];
    case 3
        load toast_demo5 prior_sum
        prm.regul.prior.refimg = [prior_sum; prior_sum];
    case 4
        load toast_demo5 prior_part
        prm.regul.prior.refimg = [prior_part; prior_part];
    case 5
        load toast_demo5 prior_rot
        prm.regul.prior.refimg = [prior_rot; prior_rot];
end
setappdata(handles.figure1,'prm',prm);
showref(handles,prm);


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


% --- Executes on slider movement.
function slider8_Callback(hObject, eventdata, handles)
% hObject    handle to slider8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
prm = getappdata(handles.figure1,'prm');
prm.regul.prior.smooth = get(hObject,'Value');
set(handles.text28,'String',num2str(prm.regul.prior.smooth));
setappdata(handles.figure1,'prm',prm);
showref(handles,prm);


% --- Executes during object creation, after setting all properties.
function slider8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider9_Callback(hObject, eventdata, handles)
% hObject    handle to slider9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
prm = getappdata(handles.figure1,'prm');
prm.regul.prior.threshold = get(hObject,'Value');
set(handles.text30,'String',num2str(prm.regul.prior.threshold));
setappdata(handles.figure1,'prm',prm);
showref(handles,prm);


% --- Executes during object creation, after setting all properties.
function slider9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end




% --- Executes on slider movement.
function slider10_Callback(hObject, eventdata, handles)
% hObject    handle to slider10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
prm = getappdata(handles.figure1,'prm');
prm.regul.tau = 10^(get(hObject,'Value'));
set(handles.text32,'String',num2str(prm.regul.tau,'%0.3e'));
setappdata(handles.figure1,'prm',prm);


% --- Executes during object creation, after setting all properties.
function slider10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on selection change in popupmenu3.
function popupmenu3_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu3
prm = getappdata(handles.figure1,'prm');

switch get(hObject,'Value')
    case 1
        set(handles.slider11,'Visible','off');
        set(handles.text34,'Visible','off');
        set(handles.text35,'Visible','off');
        prm.regul.method = 'TK1';
    case 2
        prm.regul.method = 'TV';
        set(handles.text34,'Visible','on','String','beta');
        set(handles.text35,'Visible','on','String',num2str(prm.regul.tv.beta));
        set(handles.slider11,'Visible','on','Value',log10(prm.regul.tv.beta));
    case 3
        prm.regul.method = 'Huber';
        set(handles.text34,'Visible','on','String','eps');
        set(handles.text35,'Visible','on','String',num2str(prm.regul.huber.eps));
        set(handles.slider11,'Visible','on','Value',log10(prm.regul.huber.eps));
end
setappdata(handles.figure1,'prm',prm);


% --- Executes during object creation, after setting all properties.
function popupmenu3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider11_Callback(hObject, eventdata, handles)
% hObject    handle to slider11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
prm = getappdata(handles.figure1,'prm');
v = get(hObject,'Value');
switch prm.regul.method
    case 'TV'
        prm.regul.tv.beta = 10^v;
        set(handles.text35,'String',num2str(prm.regul.tv.beta));
    case 'Huber'
        prm.regul.huber.eps = 10^v;
        set(handles.text35,'String',num2str(prm.regul.huber.eps));
end
setappdata(handles.figure1,'prm',prm);


% --- Executes during object creation, after setting all properties.
function slider11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

