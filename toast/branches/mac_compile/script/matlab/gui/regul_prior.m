function varargout = regul_prior(varargin)
%REGUL_PRIOR M-file for regul_prior.fig
%      REGUL_PRIOR, by itself, creates a new REGUL_PRIOR or raises the existing
%      singleton*.
%
%      H = REGUL_PRIOR returns the handle to a new REGUL_PRIOR or the handle to
%      the existing singleton*.
%
%      REGUL_PRIOR('Property','Value',...) creates a new REGUL_PRIOR using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to regul_prior_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      REGUL_PRIOR('CALLBACK') and REGUL_PRIOR('CALLBACK',hObject,...) call the
%      local function named CALLBACK in REGUL_PRIOR.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help regul_prior

% Last Modified by GUIDE v2.5 08-Mar-2008 00:30:25

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @regul_prior_OpeningFcn, ...
                   'gui_OutputFcn',  @regul_prior_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
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


% --- Executes just before regul_prior is made visible.
function regul_prior_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for regul_prior
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% store the parameter instance
prm = varargin{1};
setappdata(hObject,'prm',prm);
initgui(handles,prm);
updategui(handles,prm);

% UIWAIT makes regul_prior wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = regul_prior_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
%varargout{1} = handles.output;
if length(handles) > 0
    prm = getappdata(handles.figure1,'prm');
    prm = updateprm(handles,prm);
    varargout{1} = prm;
    delete(handles.figure1);
else
    varargout{1} = 0;
end

function regul_tv_beta_Callback(hObject, eventdata, handles)
% hObject    handle to regul_tv_beta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of regul_tv_beta as text
%        str2double(get(hObject,'String')) returns contents of regul_tv_beta as a double


% --- Executes during object creation, after setting all properties.
function regul_tv_beta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to regul_tv_beta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function prior_refimg_Callback(hObject, eventdata, handles)
% hObject    handle to prior_refimg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of prior_refimg as text
%        str2double(get(hObject,'String')) returns contents of prior_refimg as a double


% --- Executes during object creation, after setting all properties.
function prior_refimg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to prior_refimg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


function prm = getprm_Callback(hObject,eventdata,handles)
prm = getappdata(handles.figure1,'prm');


function initgui(handles,prm)
set(handles.uipanel1,'Title',[prm.method ' parameters']);
set(handles.text1,'Visible','off');
set(handles.regul_tv_beta,'Visible','off');
set(handles.text5,'Visible','off');
set(handles.regul_huber_eps,'Visible','off');
switch prm.method
    case 'TV'
        set(handles.text1,'Visible','on');
        set(handles.regul_tv_beta,'Visible','on');
    case 'Huber'
        set(handles.text5,'Visible','on');
        set(handles.regul_huber_eps,'Visible','on');
end


function updategui(handles,prm)
switch prm.method
    case 'TV'
        set(handles.regul_tv_beta,'String',num2str(prm.tv.beta));
    case 'Huber'
        set(handles.regul_huber_eps,'String',num2str(prm.huber.eps));
end
updateprior(handles,prm);


function updateprior(handles,prm)
p = isfield(prm,'prior');
if p
    p = isfield(prm.prior,'disable') == false || prm.prior.disable == false;
end
set(handles.checkbox1,'Value',p);
tmp={'off' 'on'};
vis = tmp{p+1};
set(handles.text2,'Visible',vis);
set(handles.prior_refimg,'Visible',vis);
set(handles.pushbutton1,'Visible',vis);
set(handles.text3,'Visible',vis);
set(handles.text4,'Visible',vis);
set(handles.prior_smooth,'Visible',vis);
set(handles.prior_th,'Visible',vis);
if p
    if isfield(prm.prior,'refname') == false
        prm.prior.refname = '';
    end
    set(handles.prior_refimg,'String',prm.prior.refname);
    if isfield(prm.prior,'smooth') == false
        prm.prior.smooth = 1;
    end
    set(handles.prior_smooth,'String',num2str(prm.prior.smooth));
    if isfield(prm.prior,'threshold') == false
        prm.prior.threshold = 0.1;
    end
    set(handles.prior_th,'String',num2str(prm.prior.threshold));
end


function p = updateprm(handles,prm)
p = prm;
switch prm.method
    case 'TV'
        p.tv.beta = str2num(get(handles.regul_tv_beta,'String'));
    case 'Huber'
        p.huber.eps = str2num(get(handles.regul_huber_eps,'String'));
end
doprior = get(handles.checkbox1,'Value');
if doprior
    p.prior.refname = get(handles.prior_refimg,'String');
    if length(p.prior.refname) == 0
        p.prior = rmfield(p.prior,'refname');
    end
    p.prior.smooth = str2num(get(handles.prior_smooth,'String'));
    p.prior.threshold = str2num(get(handles.prior_th,'String'));
    if isfield(p.prior,'disable') && p.prior.disable == true
        doprior = false;
    end
end
if doprior == false && isfield(p,'prior')
    p = rmfield(p,'prior');
end


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiresume(handles.figure1);



function prior_smooth_Callback(hObject, eventdata, handles)
% hObject    handle to prior_smooth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of prior_smooth as text
%        str2double(get(hObject,'String')) returns contents of prior_smooth as a double


% --- Executes during object creation, after setting all properties.
function prior_smooth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to prior_smooth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function prior_th_Callback(hObject, eventdata, handles)
% hObject    handle to prior_th (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of prior_th as text
%        str2double(get(hObject,'String')) returns contents of prior_th as a double


% --- Executes during object creation, after setting all properties.
function prior_th_CreateFcn(hObject, eventdata, handles)
% hObject    handle to prior_th (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function regul_huber_eps_Callback(hObject, eventdata, handles)
% hObject    handle to regul_huber_eps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of regul_huber_eps as text
%        str2double(get(hObject,'String')) returns contents of regul_huber_eps as a double


% --- Executes during object creation, after setting all properties.
function regul_huber_eps_CreateFcn(hObject, eventdata, handles)
% hObject    handle to regul_huber_eps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
prm = getappdata(handles.figure1,'prm');
enable = get(hObject,'Value');
if enable 
    if isfield(prm,'prior') == false
        prm.prior = struct ('refname','');
    end
    if isfield(prm.prior,'disable') == true
        prm.prior = rmfield (prm.prior,'disable');
    end
else
    if isfield(prm,'prior') == true
        prm.prior.disable = true;
    end
end
setappdata(handles.figure1,'prm',prm);
updateprior(handles,prm);
