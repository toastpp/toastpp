function varargout = toast_demo4(varargin)
% TOAST_DEMO4 M-file for toast_demo4.fig
%      TOAST_DEMO4, by itself, creates a new TOAST_DEMO4 or raises the existing
%      singleton*.
%
%      H = TOAST_DEMO4 returns the handle to a new TOAST_DEMO4 or the handle to
%      the existing singleton*.
%
%      TOAST_DEMO4('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TOAST_DEMO4.M with the given input arguments.
%
%      TOAST_DEMO4('Property','Value',...) creates a new TOAST_DEMO4 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before toast_demo4_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to toast_demo4_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help toast_demo4

% Last Modified by GUIDE v2.5 27-Feb-2008 13:01:42

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @toast_demo4_OpeningFcn, ...
                   'gui_OutputFcn',  @toast_demo4_OutputFcn, ...
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


% --- Executes just before toast_demo4 is made visible.
function toast_demo4_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to toast_demo4 (see VARARGIN)

% Choose default command line output for toast_demo4
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

init(handles);

% UIWAIT makes toast_demo4 wait for user response (see UIRESUME)
%uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = toast_demo4_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


%% ===========================================================
function init(handles)

prm = toastParam;
prm.fwdsolver.meshfile = 'circle25_32.msh';
prm.solver.basis.bdim = [64 64];
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
prm.initprm.mua = struct('reset','HOMOG','val',0.025);
prm.initprm.mus = struct('reset','HOMOG','val',2);
prm.initprm.ref = struct('reset','HOMOG','val',1.4);

prm.transient.callback.context = handles;
prm.transient.callback.iter = @callback_vis; % iteration callback from recon

prm.fwdsolver.hmesh = toastMesh(prm.fwdsolver.meshfile);
prm.fwdsolver.hmesh.ReadQM (prm.meas.qmfile);

prm.solver.basis.bdim = [64 64];
prm.solver.basis.hbasis = toastBasis(prm.fwdsolver.hmesh,prm.solver.basis.bdim,'Linear');

prm.regul = struct ('method','TK1', ...
     'tau',1e-5, ...
     'prior',struct ('refname','','smooth',1,'threshold',0.1));
prm.regul.basis = prm.solver.basis.hbasis;

n = prm.fwdsolver.hmesh.NodeCount;
load toast_demo1.mat
load toast_demo2.mat
prm.initprm.mua.bmua = bmua; clear bmua;
prm.initprm.mus.bmus = bmus; clear bmus;
prm.initprm.mua.mua = prm.solver.basis.hbasis.Map('B->M',prm.initprm.mua.bmua);
prm.initprm.mus.mus = prm.solver.basis.hbasis.Map('B->M',prm.initprm.mus.bmus);
prm.initprm.ref.ref = ones(n,1)*1.4;
axes(handles.axes1);
imagesc(rot90(reshape(prm.initprm.mua.bmua,prm.solver.basis.bdim)),[0.005 0.05]); axis xy equal tight off
%set(handles.axes1,'XTick',[],'XTickLabel','','YTick',[],'YTickLabel','');
axes(handles.axes2);
imagesc(rot90(reshape(prm.initprm.mus.bmus,prm.solver.basis.bdim)),[0.5 4]); axis xy equal tight off
%set(handles.axes2,'XTick',[],'XTickLabel','','YTick',[],'YTickLabel','');

setappdata(handles.figure1,'prm',prm);
disp_detectors(handles,prm);

set(handles.pushbutton1,'String','Run');
setappdata(handles.figure1,'isBusy',false);


% --- Displays the detector setup for source 0
function disp_detectors(handles,prm)

n = prm.fwdsolver.hmesh.NodeCount;
tmp = prm.solver.basis.hbasis.Map('M->B',ones(n,1));
axes(handles.axes7);
cla;imagesc(rot90(reshape(tmp,prm.solver.basis.bdim))); axis xy equal tight off
%cla;surface(reshape(tmp,prm.bx,prm.by),'EdgeColor','none'); axis tight
hold on;

qp = prm.fwdsolver.hmesh.Qpos();
qidx = 1;
x = (qp(qidx,1)/50.0*0.96+0.5)*prm.solver.basis.bdim(1)+1;
y = (qp(qidx,2)/50.0*0.96+0.5)*prm.solver.basis.bdim(2)+1;
plot3(x, y, 1, 'og', 'MarkerSize',3, 'MarkerFaceColor','green');

mp = prm.fwdsolver.hmesh.Mpos();
LL = prm.fwdsolver.hmesh.DataLinkList('matrix');
%LL = toastGetQM(prm.basis.hMesh);
idx = find(LL(:,qidx));
for i=1:length(idx)
    x = (mp(idx(i),1)/50.0*0.96+0.5)*prm.solver.basis.bdim(1)+1;
    y = (mp(idx(i),2)/50.0*0.96+0.5)*prm.solver.basis.bdim(2)+1;
    plot3(x, y, 1, 'oy', 'MarkerSize',3, 'MarkerFaceColor','yellow');
end
clear tmp
set(handles.text25,'String',['(shown is source ' num2str(qidx) ' of ' num2str(size(mp,1)) ')']);


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
function callback_vis(prm,res)
handles = prm.transient.callback.context;
axes(handles.axes3);
imagesc(rot90(reshape(res.bmua,prm.solver.basis.bdim)),[0.005 0.05]); axis xy equal tight off
%set(handles.axes3,'XTick',[],'XTickLabel','','YTick',[],'YTickLabel','');
axes(handles.axes4);
imagesc(rot90(reshape(res.bmus,prm.solver.basis.bdim)),[0.5 4]); axis xy equal tight off
%set(handles.axes4,'XTick',[],'XTickLabel','','YTick',[],'YTickLabel','');
axes(handles.axes5);
semilogy(res.of); axis tight
set(handles.axes5,'FontSize',7);
drawnow


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
        prm.meas.qmfile = 'circle25_16x16.qm';
        prm.data.lnampfile = 'fmod_ellips_16x16_100MHz.fem';
        prm.data.phasefile = 'farg_ellips_16x16_100MHz.fem';
    case 2
        prm.meas.qmfile = 'circle25_sparse.qm';
        prm.data.lnampfile = 'fmod_ellips_sparse.fem';
        prm.data.phasefile = 'farg_ellips_sparse.fem';
    case 3
        prm.meas.qmfile = 'circle25_wfan.qm';
        prm.data.lnampfile = 'fmod_ellips_wfan.fem';
        prm.data.phasefile = 'farg_ellips_wfan.fem';
    case 4
        prm.meas.qmfile = 'circle25_nfan.qm';
        prm.data.lnampfile = 'fmod_ellips_nfan.fem';
        prm.data.phasefile = 'farg_ellips_nfan.fem';
    case 5
        prm.meas.qmfile = 'circle25_wback.qm';
        prm.data.lnampfile = 'fmod_ellips_wback.fem';
        prm.data.phasefile = 'farg_ellips_wback.fem';
    case 6
        prm.meas.qmfile = 'circle25_nback.qm';
        prm.data.lnampfile = 'fmod_ellips_nback.fem';
        prm.data.phasefile = 'farg_ellips_nback.fem';
end

prm.fwdsolver.hmesh.ReadQM (prm.meas.qmfile);
setappdata(handles.figure1,'prm',prm);
disp_detectors(handles,prm);


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




