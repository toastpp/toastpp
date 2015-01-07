function varargout = recon_data_gui(varargin)
% RECON_DATA_GUI M-file for recon_data_gui.fig
%      RECON_DATA_GUI, by itself, creates a new RECON_DATA_GUI or raises the existing
%      singleton*.
%
%      H = RECON_DATA_GUI returns the handle to a new RECON_DATA_GUI or the handle to
%      the existing singleton*.
%
%      RECON_DATA_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RECON_DATA_GUI.M with the given input arguments.
%
%      RECON_DATA_GUI('Property','Value',...) creates a new RECON_DATA_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before recon_data_gui_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to recon_data_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help recon_data_gui

% Last Modified by GUIDE v2.5 03-Dec-2008 19:00:44

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @recon_data_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @recon_data_gui_OutputFcn, ...
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


% --- Executes just before recon_data_gui is made visible.
function recon_data_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to recon_data_gui (see VARARGIN)

% Choose default command line output for recon_data_gui
handles.output = hObject;

set (handles.figure1, 'name', 'Recon input data');

% Update handles structure
guidata(hObject, handles);

% note that the prm struct contains only the data subsection of the global
% parameter structure.
if nargin >= 4
    prm = varargin{1};
else
    prm = prm_setdefaults();
end
setappdata(hObject,'prm',prm);
setguifromprm(handles,prm);

% UIWAIT makes recon_data_gui wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = recon_data_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
if isstruct(handles)
    prm = getappdata(handles.figure1,'prm');
    varargout{1} = prm;
    delete(handles.figure1);
else
    varargout{1} = [];
end

function prm = prm_setdefaults()
prm = [];
prm.freq = 0;
prm.lnampfile = '';
prm.phasefile = '';
prm.useref = false;


function setguifromprm(handles,prm)
if isfield(prm,'freq')
    set(handles.edit1,'String',prm.freq);
else
    set(handles.edit1,'String','0');
end
if isfield(prm,'lnampfile')
    set(handles.edit2,'String',prm.lnampfile);
else
    set(handles.edit2,'String','');
end
if isfield(prm,'phasefile')
    set(handles.edit3,'String',prm.phasefile);
else
    set(handles.edit3,'String','');
end
if isfield(prm,'useref') && prm.useref ~= false
    haveref = true;
    if prm.useref
        set(handles.radiobutton2,'Value',1);
    else
        set(handles.radiobutton1,'Value',1);
    end
else
    haveref = false;
    set(handles.radiobutton1,'Value',1);
end
if isfield(prm,'ref')
    if isfield(prm.ref,'lnampfile')
        set(handles.edit4,'String',prm.ref.lnampfile);
    else
        set(handles.edit4,'String','');
    end
    if isfield(prm.ref,'phasefile')
        set(handles.edit5,'String',prm.ref.phasefile);
    else
        set(handles.edit5,'String','');
    end
end
showref(handles,haveref);


function showref(handles,haveref)
if haveref
    show = 'on';
    h = 32;
else
    show = 'off';
    h = 24;
end
set(handles.uipanel4,'Visible',show);
pos = get(handles.figure1,'Position');
pos(4) = h;
%set(handles.figure1,'Position',pos);





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



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
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
[fname,fpath] = uigetfile('*.fem;*.dat','Load data file (log amplitude)');
if fname
    set(handles.edit2,'String',relpath([fpath fname],pwd()));
end


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[fname,fpath] = uigetfile('*.fem;*.dat','Load data file (phase)');
if fname
    set(handles.edit3,'String',relpath([fpath fname],pwd()));
end


function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[fname,fpath] = uigetfile('*.fem;*.dat','Load data file (log amplitude)');
if fname
    set(handles.edit4,'String',relpath([fpath fname],pwd()));
end


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[fname,fpath] = uigetfile('*.fem;*.dat','Load data file (phase)');
if fname
    set(handles.edit5,'String',relpath([fpath fname],pwd()));
end


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
prm = getappdata(handles.figure1,'prm');
prm.freq = str2num(get(handles.edit1,'String'));
prm.lnampfile = get(handles.edit2,'String');
prm.phasefile = get(handles.edit3,'String');
prm.useref = logical(get(handles.radiobutton2,'Value'));
if prm.useref
    prm.ref.lnampfile = get(handles.edit4,'String');
    prm.ref.phasefile = get(handles.edit5,'String');
else
    if isfield(prm,'ref')
        rmfield(prm,'ref');
    end
end
setappdata(handles.figure1,'prm',prm);
uiresume(handles.figure1);


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(handles.figure1);


function dmode_Callback(hObject, eventdata, handles)
prm = getappdata(handles.figure1,'prm');
if hObject == handles.radiobutton1
    prm.useref = false;
else
    prm.useref = true;
end
setappdata(handles.figure1,'prm',prm);
showref(handles,prm.useref);


% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
eval('doc recon_data_gui');

