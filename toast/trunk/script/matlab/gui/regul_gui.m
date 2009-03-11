function varargout = regul_gui(varargin)
% REGUL_GUI M-file for regul_gui.fig
%      REGUL_GUI, by itself, creates a new REGUL_GUI or raises the existing
%      singleton*.
%
%      H = REGUL_GUI returns the handle to a new REGUL_GUI or the handle to
%      the existing singleton*.
%
%      REGUL_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in REGUL_GUI.M with the given input arguments.
%
%      REGUL_GUI('Property','Value',...) creates a new REGUL_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before regul_gui_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to regul_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help regul_gui

% Last Modified by GUIDE v2.5 12-Jan-2009 14:14:33

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @regul_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @regul_gui_OutputFcn, ...
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


% --- Executes just before regul_gui is made visible.
function regul_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to regul_gui (see VARARGIN)

% Choose default command line output for regul_gui
handles.output = hObject;

set(handles.figure1,'name','Regularisation parameters');

% Update handles structure
guidata(hObject, handles);

if nargin >= 4
    prm = varargin{1};
else
    prm = prm_setdefaults();
end
setappdata(hObject,'prm',prm);
setguifromprm(handles,prm);

% UIWAIT makes regul_gui wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = regul_gui_OutputFcn(hObject, eventdata, handles) 
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
prm.method = 'none';


function setguifromprm(handles,prm)
if isfield(prm,'method')
    idx = find(strcmpi({'none','TK1','TV','Huber'}, prm.method));
else
    idx = 1;
end
set (handles.popupmenu1,'Value',idx);
if isfield(prm,'tau')
    set(handles.edit1,'String',num2str(prm.tau));
end
showhide(handles);


function showhide(handles)
idx = get(handles.popupmenu1,'Value');
if idx == 1
    regvis = 'off';
else
    regvis = 'on';
end
if idx <= 2
    prmvis = 'off';
else
    prmvis = 'on';
end
reghandle = [handles.uipanel2,handles.text1,handles.edit1];
prmhandle = [handles.pushbutton4];
for i=1:length(reghandle)
    set(reghandle(i),'Visible',regvis);
end
for i=1:length(prmhandle)
    set(prmhandle(i),'Visible',prmvis);
end


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1
showhide(handles);


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


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
prm = getappdata(handles.figure1,'prm');
canclose = true;
idx = get(handles.popupmenu1,'Value');
mtd = {'none','TK1','TV','Huber'};
prm.method = mtd(idx);
prm.method = prm.method{1};
tau = str2num(get(handles.edit1,'String'));
if length(tau) > 0
    prm.tau = tau;
elseif idx > 1
    canclose = false;
end
setappdata(handles.figure1,'prm',prm);
if canclose
    uiresume(handles.figure1);
end


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(handles.figure1);


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
prm = getappdata(handles.figure1,'prm');
idx = get(handles.popupmenu1,'Value');
switch idx
    case 3
        if isfield(prm,'tv')
            tmp = regul_tv(prm.tv);
        else
            tmp = regul_tv();
        end
        if isstruct(tmp)
            prm.tv = tmp;
            setappdata(handles.figure1,'prm',prm);
        end
    case 4
        if isfield(prm,'huber')
            tmp = regul_huber(prm.huber);
        else
            tmp = regul_huber();
        end
        if isstruct(tmp)
            prm.huber = tmp;
            setappdata (handles.figure1,'prm',prm);
        end
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


