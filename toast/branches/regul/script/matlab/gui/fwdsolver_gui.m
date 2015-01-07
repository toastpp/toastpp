function varargout = fwdsolver_gui(varargin)
% FWDSOLVER_GUI M-file for fwdsolver_gui.fig
%      FWDSOLVER_GUI, by itself, creates a new FWDSOLVER_GUI or raises the existing
%      singleton*.
%
%      H = FWDSOLVER_GUI returns the handle to a new FWDSOLVER_GUI or the handle to
%      the existing singleton*.
%
%      FWDSOLVER_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FWDSOLVER_GUI.M with the given input arguments.
%
%      FWDSOLVER_GUI('Property','Value',...) creates a new FWDSOLVER_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before fwdsolver_gui_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to fwdsolver_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help fwdsolver_gui

% Last Modified by GUIDE v2.5 21-Nov-2008 19:12:57

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @fwdsolver_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @fwdsolver_gui_OutputFcn, ...
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


% --- Executes just before fwdsolver_gui is made visible.
function fwdsolver_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to fwdsolver_gui (see VARARGIN)

% Choose default command line output for fwdsolver_gui
handles.output = hObject;

set(handles.figure1,'name','Forward solver parameters');

% Update handles structure
guidata(hObject, handles);

if nargin >= 4
    prm = varargin{1};
else
    prm = prm_setdefaults();
end
setappdata(hObject,'prm',prm);
setguifromprm(handles,prm);

% UIWAIT makes fwdsolver_gui wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = fwdsolver_gui_OutputFcn(hObject, eventdata, handles) 
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
prm.meshfile = [];
prm.method = 'Direct';
prm.tol = 1e-10;


function setguifromprm(handles,prm)
if isfield(prm,'meshfile')
    set(handles.edit1,'String',prm.meshfile);
else
    set(handles.edit1,'String','');
end
if isfield(prm,'method')
    set(handles.popupmenu1,'Value',find(strcmpi(get(handles.popupmenu1,'String'),prm.method)));
end
if isfield(prm,'tol')
    set(handles.edit2,'String',num2str(prm.tol));
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


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[fname,fpath] = uigetfile('*.msh;*.opt','Load FEM mesh file');
if fname
    set(handles.edit1,'String',relpath([fpath fname],pwd()));
end


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


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


% --- Executes on button press in ok.
function ok_Callback(hObject, eventdata, handles)
% hObject    handle to ok (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
prm = getappdata(handles.figure1,'prm');
prm.meshfile = get(handles.edit1,'String');
tmp = get(handles.popupmenu1,'String');
prm.method = tmp{get(handles.popupmenu1,'Value')};
prm.tol = str2num(get(handles.edit2,'String'));
setappdata(handles.figure1,'prm',prm);
uiresume(handles.figure1);


% --- Executes on button press in cancel.
function cancel_Callback(hObject, eventdata, handles)
% hObject    handle to cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(handles.figure1);

