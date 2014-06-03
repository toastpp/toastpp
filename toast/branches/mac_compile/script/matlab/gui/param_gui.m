function varargout = param_gui(varargin)
% PARAM_GUI M-file for param_gui.fig
%      PARAM_GUI, by itself, creates a new PARAM_GUI or raises the existing
%      singleton*.
%
%      H = PARAM_GUI returns the handle to a new PARAM_GUI or the handle to
%      the existing singleton*.
%
%      PARAM_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PARAM_GUI.M with the given input arguments.
%
%      PARAM_GUI('Property','Value',...) creates a new PARAM_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before param_gui_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to param_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help param_gui

% Last Modified by GUIDE v2.5 03-Dec-2008 11:37:38

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @param_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @param_gui_OutputFcn, ...
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


% --- Executes just before param_gui is made visible.
function param_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to param_gui (see VARARGIN)

% Choose default command line output for param_gui
handles.output = hObject;

set(handles.figure1,'name','Optical parameters');

% Update handles structure
guidata(hObject, handles);

if nargin >= 4
    prm = varargin{1};
else
    prm = prm_setdefaults();
end
setappdata(hObject,'prm',prm);
setguifromprm(handles,prm);

% UIWAIT makes param_gui wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = param_gui_OutputFcn(hObject, eventdata, handles) 
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
prm.mua.reset = 'Mesh';
prm.mus.reset = 'Mesh';
prm.ref.reset = 'Mesh';


function setguifromprm(handles,prm)
if isfield(prm,'mua')
    h = [handles.popupmenu1, handles.edit1, handles.text1, handles.pushbutton3];
    updategui_prm(prm.mua,h);
    h = [handles.popupmenu2, handles.edit2, handles.text2, handles.pushbutton4];
    updategui_prm(prm.mus,h);
    h = [handles.popupmenu3, handles.edit3, handles.text3, handles.pushbutton5];
    updategui_prm(prm.ref,h);
end


    
function updategui_prm(prm,h)
idx = find(strcmpi({'HOMOG' 'NIM' 'NIM_LAST' 'MESH'}, prm.reset));
set(h(1),'Value',idx);
switch idx
    case 1
        set(h(2),'Visible','on');
        if isfield(prm,'val')
            set(h(2),'String',prm.val);
        else
            set(h(2),'String','0');
        end
        set(h(3),'String','Value [1/mm]:');
        set(h(4),'Visible','off');
    case 2
        set(h(2),'Visible','on');
        if isfield(prm,'nim')
            set(h(2),'String',prm.nim);
        else
            set(h(2),'String','');
        end
        set(h(3),'String','NIM file:');
        set(h(4),'Visible','on');
    case 3
        set(h(2),'Visible','on');
        if isfield(prm,'nim')
            set(h(2),'String',prm.nim);
        else
            set(h(2),'String','');
        end
        set(h(3),'String','NIM file:');
        set(h(4),'Visible','on');
    case 4
        set(h(2),'Visible','off');
        set(h(3),'String','');
        set(h(4),'Visible','off');
end


function p = updateprm_gui(h)
opt = {'HOMOG' 'NIM' 'NIM_LAST' 'MESH'};
idx = get(h(1),'Value');
p.reset = opt{idx};
switch idx
    case 1
        p.val = str2num(get(h(2),'String'));
    case 2
        p.nim = get(h(2),'String');
    case 3
        p.nim = get(h(2),'String');
end


% --- Executes on button press in ok.
function ok_Callback(hObject, eventdata, handles)
% hObject    handle to ok (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
prm = getappdata(handles.figure1,'prm');
prm.mua = updateprm_gui([handles.popupmenu1,handles.edit1]);
prm.mus = updateprm_gui([handles.popupmenu2,handles.edit2]);
prm.ref = updateprm_gui([handles.popupmenu3,handles.edit3]);
setappdata(handles.figure1,'prm',prm);
uiresume(handles.figure1);


% --- Executes on button press in cancel.
function cancel_Callback(hObject, eventdata, handles)
% hObject    handle to cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(handles.figure1);


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1
prm = getappdata(handles.figure1,'prm');
opt = {'HOMOG' 'NIM' 'NIM_LAST' 'MESH'};
prm.mua.reset = opt{get(hObject,'Value')};
setappdata(handles.figure1,'prm',prm);
updategui_prm(prm.mua,[handles.popupmenu1, handles.edit1, handles.text1, handles.pushbutton3]);



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


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[fname,fpath] = uigetfile('*.nim','Load NIM image file');
if fname
    set(handles.edit1,'String',relpath([fpath fname],pwd()));
end


% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2
prm = getappdata(handles.figure1,'prm');
opt = {'HOMOG' 'NIM' 'NIM_LAST' 'MESH'};
prm.mus.reset = opt{get(hObject,'Value')};
setappdata(handles.figure1,'prm',prm);
updategui_prm(prm.mus,[handles.popupmenu2, handles.edit2, handles.text2, handles.pushbutton4]);


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


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[fname,fpath] = uigetfile('*.nim','Load NIM image file');
if fname
    set(handles.edit2,'String',relpath([fpath fname],pwd()));
end


% --- Executes on selection change in popupmenu3.
function popupmenu3_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu3
prm = getappdata(handles.figure1,'prm');
opt = {'HOMOG' 'NIM' 'NIM_LAST' 'MESH'};
prm.ref.reset = opt{get(hObject,'Value')};
setappdata(handles.figure1,'prm',prm);
updategui_prm(prm.ref,[handles.popupmenu3, handles.edit3, handles.text3, handles.pushbutton5]);


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


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[fname,fpath] = uigetfile('*.nim','Load NIM image file');
if fname
    set(handles.edit3,'String',relpath([fpath fname],pwd()));
end

