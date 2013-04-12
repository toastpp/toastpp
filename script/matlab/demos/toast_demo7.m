function varargout = toast_demo7(varargin)
% TOAST_demo7 M-file for toast_demo7.fig
%      TOAST_demo7, by itself, creates a new TOAST_demo7 or raises the existing
%      singleton*.
%
%      H = TOAST_demo7 returns the handle to a new TOAST_demo7 or the handle to
%      the existing singleton*.
%
%      TOAST_demo7('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TOAST_demo7.M with the given input arguments.
%
%      TOAST_demo7('Property','Value',...) creates a new TOAST_demo7 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before toast_demo7_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to toast_demo7_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help toast_demo7

% Last Modified by GUIDE v2.5 08-Apr-2008 12:49:40

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @toast_demo7_OpeningFcn, ...
                   'gui_OutputFcn',  @toast_demo7_OutputFcn, ...
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


% --- Executes just before toast_demo7 is made visible.
function toast_demo7_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to toast_demo7 (see VARARGIN)

% Choose default command line output for toast_demo7
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

init(handles);

% UIWAIT makes toast_demo7 wait for user response (see UIRESUME)
%uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = toast_demo7_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% ===========================================================
function init(handles)

rand('state',14);
prm = toastParam;
prm.user.maxch = 4;
prm.user.maxlambda = 4;
prm.user.nch = 2;
prm.user.nlambda = 2;
prm.user.extinct = [[0.01 0.02 0.03 0.035]; [0.035 0.015 0.01 0.02]; [0.03 0.02 0.02 0.03]; [0.015 0.03 0.035 0.025]]; % extinction coefficient for both chromophores at both wavelengths
prm.solver.basis.bdim = [256 256]; %[64 64];
prm.fwdsolver.meshfile = 'circle25_32.msh';
prm.meas.qmfile = 'circle25_16x16.qm';
prm.meas.src = struct('freq',100,'type','Neumann','prof','Gaussian','width',2);
prm.meas.det = struct('prof','Gaussian','width',2);
prm.fwdsolver.method = 'Direct';
prm.solver.method = 'LM';
prm.solver.tol = 1e-8;
prm.solver.step0 = 100;
prm.solver.krylov = struct('method','gmres','tol',1e-2,'maxit',100);
prm.regul.method = 'None';
prm.initprm.mus = struct('reset','HOMOG','val',2);
prm.initprm.ref = struct('reset','HOMOG','val',1.4);

prm.user.edit_extinct.ch = 1;
prm.user.edit_extinct.lambda = 1;

prm.transient.callback.context = handles;
prm.transient.callback.iter = @callback_vis; % iteration callback from recon

prm.fwdsolver.hmesh = toastMesh(prm.fwdsolver.meshfile);
prm.fwdsolver.hmesh.ReadQM (prm.meas.qmfile);

prm.solver.basis.hbasis = toastBasis(prm.fwdsolver.hmesh,prm.solver.basis.bdim,'Linear');
n = prm.fwdsolver.hmesh.NodeCount;
%load toast_demo1.mat
%load toast_demo2.mat
%prm.initprm.mua.bmua = bmua; clear bmua;
%prm.initprm.mus.bmus = bmus; clear bmus;
%prm.initprm.mua.mua = prm.solver.basis.hbasis.Map('B->M',prm.initprm.mua.bmua);
%prm.initprm.mus.mus = prm.solver.basis.hbasis.Map('B->M',prm.initprm.mus.bmus);
%prm.initprm.ref.ref = ones(n,1)*1.4;
prm.initprm.mus.mus = ones(n,1)*2;
prm.initprm.ref.ref = ones(n,1)*1.4;
prm.initprm.mus.bmus = prm.solver.basis.hbasis.Map('M->B',prm.initprm.mus.mus);


prm.user.hTgt = zeros(prm.user.maxch,1);
prm.user.hMuaTgt = zeros(prm.user.maxlambda,1);
setappdata(handles.figure1,'prm',prm);
rolltgt (handles);
updategui (handles);

set(handles.pushbutton1,'String','Run');
setappdata(handles.figure1,'isBusy',false);


% ===========================================================
function show_C_images (handles)
prm = getappdata(handles.figure1,'prm');
for i=1:4
    if i <= prm.user.nch
        vis = 'on';
    else
        vis = 'off';
    end
    eval (['set(handles.text' num2str(i) ',''Visible'',''' vis ''')']);
    eval (['set(handles.text' num2str(i+4) ',''Visible'',''' vis ''')']);
    eval (['set(handles.text' num2str(i+12) ',''Visible'',''' vis ''')']);

    if i <= prm.user.nch
        eval(['axes(handles.axes' num2str(i) ')']);
        h = imagesc(reshape(prm.user.bC_tgt(i,:), prm.solver.basis.bdim(1), ...
                            prm.solver.basis.bdim(2)),[0,prm.user.bCmax]);
        prm.user.hTgt(i) = h;
    elseif prm.user.hTgt(i) ~= 0
        delete(prm.user.hTgt(i));
        prm.user.hTgt(i) = 0;
    end
    eval (['set(handles.axes' num2str(i) ',''Visible'',''off'')']);
    eval (['set(handles.axes' num2str(i+4) ',''Visible'',''off'')']);
    setappdata(handles.figure1,'prm',prm);
end


% ===========================================================
function show_mua_images (handles)
prm = getappdata(handles.figure1,'prm');
for i=1:prm.user.nlambda
    img(i,:) = zeros(1,prod(prm.solver.basis.bdim));
    for j=1:prm.user.nch
        img(i,:) = img(i,:) + prm.user.bC_tgt(j,:) * prm.user.extinct(j,i);
    end
end
muamax = max(max(img));
for i=1:4
    if i <= prm.user.nlambda
        vis = 'on';
    else
        vis = 'off';
    end
    eval (['set(handles.text' num2str(i+8) ',''Visible'',''' vis ''')']);
    if i <= prm.user.nlambda
        eval(['axes(handles.axes' num2str(i+8) ')']);
        h = imagesc(reshape(img(i,:), prm.solver.basis.bdim(1), ...
                            prm.solver.basis.bdim(2)),[0,muamax]);
        prm.user.hMuaTgt(i) = h;
    elseif prm.user.hMuaTgt(i) ~= 0
        delete(prm.user.hMuaTgt(i));
        prm.user.hMuaTgt(i) = 0;
    end
    eval (['set(handles.axes' num2str(i+8) ',''Visible'',''off'')']);
    setappdata(handles.figure1,'prm',prm);
end
clear img;


% ===========================================================
function showExtinct (handles)
prm = getappdata(handles.figure1,'prm');
col = ['b' 'r' 'g' 'k'];
axes(handles.axes13);
h = plot(prm.user.extinct(1,1:prm.user.nlambda),['-o' col(1)]);
set(h,'ButtonDownFcn','toast_demo7(''axes13_ButtonDownFcn'',gcbo,1,guidata(gcbo))');
hold on
for i=2:prm.user.nch
    h = plot(prm.user.extinct(i,1:prm.user.nlambda),['-o' col(i)]);
    set(h,'ButtonDownFcn',['toast_demo7(''axes13_ButtonDownFcn'',gcbo,' num2str(i) ',guidata(gcbo))']);
end
hold off


% ===========================================================
function updategui (handles)
prm = getappdata(handles.figure1,'prm');
set(handles.popupmenu1,'Value',prm.user.nch);
set(handles.popupmenu2,'Value',prm.user.nlambda);
set(handles.edit1,'String',num2str(prm.user.extinct(prm.user.edit_extinct.ch,prm.user.edit_extinct.lambda)));
show_C_images(handles);
show_mua_images(handles);
showExtinct(handles);


% ===========================================================
function rolltgt (handles)
prm = getappdata(handles.figure1,'prm');
grid = prm.solver.basis.bdim;
nch = 4;  % generate targets for max. number of chromophores
blobtgt = str2num(get(handles.edit2,'String'));
for c=1:nch
    img = ones(grid);
    nblob = max(1,3+round(rand()*(blobtgt-3)));
    for b=1:nblob
        rad(b) = (0.05+rand()*0.15)*grid(1);
        count = 0;
        cont = true;
        while cont && count < 5*nblob
            count = count+1;
            cont = false;
            cnt(b,:) = rand(1,2).*grid;
            r = norm(cnt(b,:)-grid/2);
            if r+rad(b) >= grid(1)/2, cont=true; end
            for bb=1:b-1
                if norm(cnt(b,:)-cnt(bb,:)) < rad(b)+rad(bb), cont=true; end
            end
        end
        mag(b) = rand()*2+1;
    end
    for j=1:grid(2)
        for i=1:grid(1)
            p = [i j];
            for b=1:nblob
                if norm(p-cnt(b,:)) <= rad(b)
                    img(i,j) = img(i,j) + mag(b);
                end
            end
        end
    end
    img = reshape(img,[],1);
    bC_tgt(c,:) = img * 0.3;
    %bC_tgt(c,:) = zeros(prod(grid),1);
    %bC_tgt(c,solmask) = img(solmask) * 0.3;
end
prm.user.bC_tgt = bC_tgt;
for i=1:size(bC_tgt,1)
    tmp(i,:) = prm.solver.basis.hbasis.Map('B->S', bC_tgt(i,:));
end
prm.user.bCmin = min(min(tmp));
prm.user.bCmax = max(max(tmp));
clear tmp;

%prm.user.bCmin = min(min(bC_tgt(:,solmask)));
%prm.user.bCmax = max(max(bC_tgt(:,solmask)));
setappdata(handles.figure1,'prm',prm);


% ===========================================================
% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if getappdata(handles.figure1,'isBusy') == false
    prm = getappdata(handles.figure1,'prm');
    setappdata(handles.figure1,'isBusy',true);
    set(handles.pushbutton1,'String','Stop (may be delayed)');
    % generate target data at all wavelengths
    prm.data = genfwddata (prm, prm.user.extinct, prm.user.bC_tgt);
    % and reconstruct concentrations from homogeneous intial guesses
    toastReconMultispectralCW(prm);
    setappdata(handles.figure1,'isBusy',false);
    set(handles.pushbutton1,'String','Run');
else
    global ITR_TERM
    ITR_TERM = true;  % send termination flag
end


% ===========================================================
% Display reconstruction results for current iteration
function callback_vis(handles,res)
prm = getappdata(handles.figure1,'prm');
for i=1:prm.user.nch
    eval(['axes(handles.axes' num2str(i+4) ')']);
    h = imagesc(reshape(res.bC(i,:), prm.solver.basis.bdim),[0,prm.user.bCmax]);
    eval (['set(handles.axes' num2str(i+4) ',''Visible'',''off'')']);
end
axes(handles.axes14);
semilogy([0:max(1,length(res.of)-1)],res.of); axis tight
set(handles.axes14,'FontSize',7);
drawnow


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1
prm = getappdata(handles.figure1,'prm');
prm.user.nch = get(hObject,'Value');
setappdata(handles.figure1,'prm',prm);
updategui(handles);


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
prm.user.nlambda = get(hObject,'Value');
setappdata(handles.figure1,'prm',prm);
updategui(handles);


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




% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
rolltgt (handles);
updategui(handles);


function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double
prm = getappdata(handles.figure1,'prm');
e = str2num(get(hObject,'String'));
prm.user.extinct(prm.user.edit_extinct.ch,prm.user.edit_extinct.lambda) = e;
setappdata(handles.figure1,'prm',prm);
updategui(handles);


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


function axes13_ButtonDownFcn(hObject, eventdata, handles)
prm = getappdata(handles.figure1,'prm');
mouse = get(gca,'currentpoint');
x = mouse(1,1);
ch = eventdata;
lambda = round(x);
prm.user.edit_extinct.ch = ch;
prm.user.edit_extinct.lambda = lambda;
set(handles.text34,'String',sprintf('Chromophore %d, wavelength %0.0f',ch,lambda));
set(handles.edit1,'String',num2str(prm.user.extinct(ch,lambda)));
setappdata(handles.figure1,'prm',prm);



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


function mua = GetMua (prm, extinct, C)
mua = zeros(size(C,2),1);
for i = 1:prm.user.nch
    mua = mua + C(i,:)'*extinct(i);
end


function data = genfwddata(prm, extinct, bC)
qvec = real(prm.fwdsolver.hmesh.Qvec(prm.meas.src.type, prm.meas.src.prof, prm.meas.src.width));
mvec = real(prm.fwdsolver.hmesh.Mvec(prm.meas.det.prof, prm.meas.det.width));
data = [];
for i=1:prm.user.nlambda
    bmua = GetMua(prm,extinct(:,i), bC);
    mua = prm.solver.basis.hbasis.Map ('G->M', bmua);
    proj = toastProject (prm.fwdsolver.hmesh, mua, prm.initprm.mus.mus, ...
        prm.initprm.ref.ref, 0, qvec, mvec, 'Direct', 1e-10);
    data = [data; proj];
end