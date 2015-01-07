function varargout = regul_tv(varargin)
% REGUL_TV M-file for regul_tv.fig
%      REGUL_TV, by itself, creates a new REGUL_TV or raises the existing
%      singleton*.
%
%      H = REGUL_TV returns the handle to a new REGUL_TV or the handle to
%      the existing singleton*.
%
%      REGUL_TV('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in REGUL_TV.M with the given input arguments.
%
%      REGUL_TV('Property','Value',...) creates a new REGUL_TV or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before regul_tv_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to regul_tv_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help regul_tv

% Last Modified by GUIDE v2.5 12-Jan-2009 15:16:34

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @regul_tv_OpeningFcn, ...
                   'gui_OutputFcn',  @regul_tv_OutputFcn, ...
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


% --- Executes just before regul_tv is made visible.
function regul_tv_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to regul_tv (see VARARGIN)

% Choose default command line output for regul_tv
handles.output = hObject;

set(handles.figure1,'name','TV regularisation parameters');

% Update handles structure
guidata(hObject, handles);

if (nargin >= 4)
    prm = varargin{1};
else
    prm = prm_setdefaults();
end
setappdata(hObject,'prm',prm);
setguifromprm(handles,prm);

% UIWAIT makes regul_tv wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = regul_tv_OutputFcn(hObject, eventdata, handles) 
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
prm.beta = 1e-2;


function setguifromprm(handles,prm)
if isfield(prm,'beta')
    set(handles.edit1, 'String', num2str(prm.beta));
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
prm = getappdata(handles.figure1,'prm');
canclose = true;
beta = str2num (get(handles.edit1,'String'));
if length(beta) > 0
    prm.beta = beta;
else
    canclose = false;
end
setappdata (handles.figure1,'prm',prm);
if canclose
    uiresume (handles.figure1);
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


