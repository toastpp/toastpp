function varargout = toast_demo1a(varargin)
% TOAST_demo1a M-file for toast_demo1a.fig
%      TOAST_demo1a, by itself, creates a new TOAST_demo1a or raises the existing
%      singleton*.
%
%      H = TOAST_demo1a returns the handle to a new TOAST_demo1a or the handle to
%      the existing singleton*.
%
%      TOAST_demo1a('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TOAST_demo1a.M with the given input arguments.
%
%      TOAST_demo1a('Property','Value',...) creates a new TOAST_demo1a or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before toast_demo1a_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to toast_demo1a_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help toast_demo1a

% Last Modified by GUIDE v2.5 24-Apr-2008 17:50:19

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @toast_demo1a_OpeningFcn, ...
                   'gui_OutputFcn',  @toast_demo1a_OutputFcn, ...
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


% --- Executes just before toast_demo1a is made visible.
function toast_demo1a_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to toast_demo1a (see VARARGIN)

% Choose default command line output for toast_demo1a
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

init(handles);

% UIWAIT makes toast_demo1a wait for user response (see UIRESUME)
%uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = toast_demo1a_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


%% ===========================================================
function init(handles)

prm.basis.bdim = [64 64];
prm.blen = prod(prm.basis.bdim);
prm.basis.hMesh = toastMesh('circle25_32.msh');
prm.basis.hMesh.ReadQM('circle25_1x32.qm');
prm.basis.hBasis = toastBasis(prm.basis.hMesh,prm.basis.bdim,'Linear');
prm.qvec = real(prm.basis.hMesh.Qvec('Isotropic','Gaussian',2));
prm.mvec = real(prm.basis.hMesh.Mvec('Gaussian',2,1.4));
prm.mvec = prm.mvec(:,[4,8,12,16]);   % only run with 4 detectors
%prm.smask = toastSolutionMask(prm.basis.hBasis);
n = prm.basis.hMesh.NodeCount();
prm.mua = ones(n,1)*0.01;
prm.mus = ones(n,1)*1;
prm.ref = ones(n,1)*1.4;
prm.bkg = 0;
prm.meas.src.freq = 0;
prm.bmua = prm.basis.hBasis.Map('M->B',prm.mua);
prm.bmus = prm.basis.hBasis.Map('M->B',prm.mus);

prm.trange = 5000;
prm.nstep = 100;
prm.theta = 0.75;

setappdata(handles.figure1,'prm',prm);
showprm(handles)
showqm(handles)

% ===========================================================
% Display parameter images

function showprm (handles)

prm = getappdata(handles.figure1,'prm');

axes(handles.axes1);
imagesc(reshape(prm.bmua,prm.basis.bdim(1),prm.basis.bdim(2)),[0 0.02]);
axis equal tight;
colorbar('location','SouthOutside','position',[0.0439 0.08 0.3672 0.05],'FontSize',7);
set(gca,'Position',[2.648 2.755 30.2 11.615],'Visible','off');

axes(handles.axes2);
imagesc(reshape(prm.bmus,prm.basis.bdim(1),prm.basis.bdim(2)),[0 2]);
axis equal tight;
colorbar('location','SouthOutside','position',[0.4408 0.08 0.3672 0.05],'FontSize',7);
set(gca,'Position',[33.6 2.755 30.2 11.615],'Visible','off');


% ===========================================================

function showqm (handles)

prm = getappdata(handles.figure1,'prm');
n = prm.basis.hMesh.NodeCount();
tmp = prm.basis.hBasis.Map('M->B',ones(n,1));

midx = [4 8 12 16];
M = prm.basis.hMesh.Mpos();
M = M(midx,:);   % only use 4 detectors
Q = prm.basis.hMesh.Qpos();
axes(handles.axes7);
cla; surface(reshape(tmp,prm.basis.bdim(1),prm.basis.bdim(2)),'EdgeColor','none');

hold on
p = (Q(1,:)./55 + 0.5) .* prm.basis.bdim;
plot3 (p(1),p(2),1,'ow','MarkerSize',7,'MarkerFaceColor','none');
p = (Q(1,:)./70 + 0.5) .* prm.basis.bdim;
text(p(1)-prm.basis.bdim(1)/20,p(2),1,'Q','Color','white');

mk = ['b' 'r' 'g' 'm'];
col = {'blue' 'red' 'green' 'magenta'};
for i=1:size(M,1)
    p = (M(i,:)./55 + 0.5) .* prm.basis.bdim;
    plot3 (p(1),p(2),1,['o' mk(i)],'MarkerSize',7, 'MarkerFaceColor',col{i});
    p = (M(i,:)./70 + 0.5) .* prm.basis.bdim;
    text(p(1)-prm.basis.bdim(1)/20,p(2),1,['D' num2str(i)],'Color','white');
end
axis equal tight;


% ===========================================================
function run_tpsf (handles)

prm = getappdata(handles.figure1,'prm');

theta = prm.theta;
nstep = prm.nstep;
dt = prm.trange/nstep;

% Set up the required FEM matrices
smat = real(dotSysmat(prm.basis.hMesh,prm.mua,prm.mus,prm.ref,0));
mmat = prm.basis.hMesh.Massmat();
K0 = -(smat * (1-theta) - mmat * 1/dt);            % backward difference matrix
K1 = smat * theta + mmat * 1/dt;                   % forward difference matrix
mvecT = prm.mvec.';

% initial condition
phi = prm.qvec;

% loop over time steps
[L U] = lu(K1);
t = [1:nstep] * dt;
ndisp = max(1,round(nstep/50));
j = ndisp-1;

tic;
for i=1:nstep
    q = K0 * phi;
    phi = U\(L\q);
    phi = max(phi,0);
    gamma(i,:) = mvecT * phi;
    gamma(i,:) = max (gamma(i,:),1e-20);
    lgamma(i,:) = log(gamma(i,:));

    j = j+1;
    if j == ndisp || i == nstep
        axes(handles.axes5);
        img = prm.basis.hBasis.Map ('M->G',log(max(phi,1e-25)));
        img = reshape(img, prm.basis.bdim(1), prm.basis.bdim(2))';
        imagesc(img,[-25 -5]); axis equal tight;
        set (gca,'Visible','off');
        axes(handles.axes3);
        plot(t(1:i),gamma(:,1)./max(gamma(:,1)));
        hold on
        plot(t(1:i),gamma(:,2)./max(gamma(:,2)),'r');
        plot(t(1:i),gamma(:,3)./max(gamma(:,3)),'g');
        plot(t(1:i),gamma(:,4)./max(gamma(:,4)),'m');
        axis tight;
        set(gca,'XTick',[],'XTickLabel','');
        hold off
    
        axes(handles.axes4);
        plot(t(1:i),lgamma(:,1));
        hold on
        plot(t(1:i),lgamma(:,2),'r');
        plot(t(1:i),lgamma(:,3),'g');
        plot(t(1:i),lgamma(:,4),'m');
        axis tight;
        hold off
        drawnow
        j = 0;
    end
end
toc

% ===========================================================
function setbackground(handles,bkg)
prm = getappdata(handles.figure1,'prm');
if bkg ~= prm.bkg
    prm.bkg = bkg;
    
    % update parameter background
    switch bkg
        case 0 % flat
            n = prm.basis.hMesh.NodeCount();
            prm.mua = ones(n,1)*0.01;
            prm.mus = ones(n,1)*1;
            prm.bmua = prm.basis.hBasis.Map('M->B',prm.mua);
            prm.bmus = prm.basis.hBasis.Map('M->B',prm.mus);
        case 1 % 'blobs'
            load toast_demo1.mat
            prm.bmua = bmua * 0.01/0.025; clear bmua;
            prm.bmus = bmus * 1/2; clear bmus;
            prm.mua = prm.basis.hBasis.Map('B->M',prm.bmua);
            prm.mus = prm.basis.hBasis.Map('B->M',prm.bmus);
    end
    setappdata(handles.figure1,'prm',prm);
    showprm(handles);
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
run_tpsf (handles);


% --- Executes on button press in radiobutton1.
function radiobutton1_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton1
prm = getappdata(handles.figure1,'prm');
if prm.bkg ~= 0
    set(handles.radiobutton2,'Value',0);
    setbackground(handles,0);
end


% --- Executes on button press in radiobutton2.
function radiobutton2_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton2
prm = getappdata(handles.figure1,'prm');
if prm.bkg ~= 1
    set(handles.radiobutton1,'Value',0);
    setbackground(handles,1);
end



% --- Executes on slider movement.
function slider8_Callback(hObject, eventdata, handles)
% hObject    handle to slider8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
prm = getappdata(handles.figure1,'prm');
prm.trange = get(hObject,'Value')*1e3;
set(handles.text33,'String',num2str(prm.trange));
setappdata(handles.figure1,'prm',prm);


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
prm.nstep = get(hObject,'Value');
set(handles.text35,'String',num2str(round(prm.nstep)));
setappdata(handles.figure1,'prm',prm);


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
prm.theta = get(hObject,'Value');
set(handles.text37,'String',num2str(prm.theta));
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


