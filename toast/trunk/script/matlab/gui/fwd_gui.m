function varargout = fwd_gui(varargin)
% FWD_GUI M-file for fwd_gui.fig
%      FWD_GUI, by itself, creates a new FWD_GUI or raises the existing
%      singleton*.
%
%      H = FWD_GUI returns the handle to a new FWD_GUI or the handle to
%      the existing singleton*.
%
%      FWD_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FWD_GUI.M with the given input arguments.
%
%      FWD_GUI('Property','Value',...) creates a new FWD_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before fwd_gui_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to fwd_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help fwd_gui

% Last Modified by GUIDE v2.5 19-Feb-2008 12:44:58

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @fwd_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @fwd_gui_OutputFcn, ...
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


% --- Executes just before fwd_gui is made visible.
function fwd_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to fwd_gui (see VARARGIN)

% Choose default command line output for fwd_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

prm = toastParam;
if nargin >= 4
    prmfile = varargin{1};
    set(handles.prmfile,'String',prmfile);
    prm = toastParam(prmfile);
    updategui(handles,prm);
    setappdata(hObject,'prmfile',prmfile);
end
setappdata(hObject,'prm',prm);

% UIWAIT makes fwd_gui wait for user response (see UIRESUME)
%uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = fwd_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function meshname_Callback(hObject, eventdata, handles)
% hObject    handle to meshname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of meshname as text
%        str2double(get(hObject,'String')) returns contents of meshname as a double


% --- Executes during object creation, after setting all properties.
function meshname_CreateFcn(hObject, eventdata, handles)
% hObject    handle to meshname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function qmname_Callback(hObject, eventdata, handles)
% hObject    handle to qmname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of qmname as text
%        str2double(get(hObject,'String')) returns contents of qmname as a double


% --- Executes during object creation, after setting all properties.
function qmname_CreateFcn(hObject, eventdata, handles)
% hObject    handle to qmname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function prmfile_Callback(hObject, eventdata, handles)
% hObject    handle to prmfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of prmfile as text
%        str2double(get(hObject,'String')) returns contents of prmfile as a double


% --- Executes during object creation, after setting all properties.
%function prmfile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to prmfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
%if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%    set(hObject,'BackgroundColor','white');
%end


% --- Update GUI elements from parameters
function updategui(handles,prm)
set(handles.meshname,'String',prm.fwdsolver.meshfile);
set(handles.qmname,'String',prm.meas.qmfile);
set(handles.lnampname,'String',prm.data.lnampfile);
set(handles.phasename,'String',prm.data.phasefile);
set(handles.freq,'String',prm.data.freq);
set(handles.srctype,'Value',find(strcmpi(get(handles.srctype,'String'),prm.meas.src.type)));
set(handles.srcprof,'Value',find(strcmpi(get(handles.srcprof,'String'),prm.meas.src.prof)));
set(handles.srcwidth,'String',prm.meas.src.width);
set(handles.mprof,'Value',find(strcmpi(get(handles.mprof,'String'),prm.meas.det.prof)));
set(handles.mwidth,'String',prm.meas.det.width);
idx = find(strcmpi(get(handles.linsolver,'String'),prm.fwdsolver.method));
if length(idx) == 0
    idx = 1;
end
set(handles.linsolver,'Value',idx);
set(handles.lintol,'String',prm.fwdsolver.tol);
updategui_initprm(prm.initprm.mua, handles.prm_mua, handles.val_mua);
updategui_initprm(prm.initprm.mus, handles.prm_mus, handles.val_mus);
updategui_initprm(prm.initprm.ref, handles.prm_ref, handles.val_ref);


% --- Update parameters from GUI elements
function prm = updateprm(handles)
prm = getappdata(handles.figure1,'prm');
prm.fwdsolver.meshfile = get(handles.meshname,'String');
prm.data.lnampfile = get(handles.lnampname,'String');
prm.data.phasefile = get(handles.phasename,'String');
prm.data.freq = str2num(get(handles.freq,'String'));
prm.meas.qmfile = get(handles.qmname,'String');
tmp = get(handles.srctype,'String');
prm.meas.src.type = tmp{get(handles.srctype,'Value')};
tmp = get(handles.srcprof,'String');
prm.meas.src.prof = tmp{get(handles.srcprof,'Value')};
prm.meas.src.width = str2num(get(handles.srcwidth,'String'));
tmp = get(handles.mprof,'String');
prm.meas.det.prof = tmp{get(handles.mprof,'Value')};
prm.meas.det.width = str2num(get(handles.mwidth,'String'));
tmp = get(handles.linsolver,'String');
prm.fwdsolver.method = tmp{get(handles.linsolver,'Value')};
prm.fwdsolver.tol = str2num(get(handles.lintol,'String'));
prm.initprm.mua = updateprm_initprm(handles.prm_mua,handles.val_mua);
prm.initprm.mus = updateprm_initprm(handles.prm_mus,handles.val_mus);
prm.initprm.ref = updateprm_initprm(handles.prm_ref,handles.val_ref);
setappdata(handles.figure1,'prm',prm);


function updategui_initprm(prm, hmethod, hval)
idx = find(strcmpi({'HOMOG' 'NIM' 'NIM_LAST' 'MESH'}, prm.reset));
set(hmethod,'Value',idx);
switch idx
    case 1
        set(hval,'String',prm.val);
    case 2
        set(hval,'String',prm.nim);
    case 3
        set(hval,'String',prm.nim);
end


function prm = updateprm_initprm(hmethod,hval)
tmp={'HOMOG' 'NIM' 'NIM_LAST' 'MESH'};
idx = get(hmethod,'Value');
prm.reset = tmp{idx};
switch idx
    case 1
        prm.val = str2num(get(hval,'String'));
    case 2
        prm.nim = get(hval,'String');
    case 3
        prm.nim = get(hval,'String');
end


function showdata(prm,mdata,pdata)
figure(1);
mesh=toastMesh(prm.fwdsolver.meshfile);
mesh.ReadQM(prm.meas.qmfile);
qpos = mesh.Qpos;
mpos = mesh.Mpos;
nq = size(qpos,1);
nm = size(mpos,1);
subplot(1,2,1); imagesc(reshape(mdata,nm,nq)); axis equal tight;
title('lnamp'); xlabel('source #'); ylabel('detector #');
subplot(1,2,2); imagesc(reshape(pdata,nm,nq)); axis equal tight;
title('phase'); xlabel('source #'); ylabel('detector #');


% --- Executes on selection change in srctype.
function srctype_Callback(hObject, eventdata, handles)
% hObject    handle to srctype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns srctype contents as cell array
%        contents{get(hObject,'Value')} returns selected item from srctype


% --- Executes during object creation, after setting all properties.
function srctype_CreateFcn(hObject, eventdata, handles)
% hObject    handle to srctype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in srcprof.
function srcprof_Callback(hObject, eventdata, handles)
% hObject    handle to srcprof (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns srcprof contents as cell array
%        contents{get(hObject,'Value')} returns selected item from srcprof


% --- Executes during object creation, after setting all properties.
function srcprof_CreateFcn(hObject, eventdata, handles)
% hObject    handle to srcprof (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function srcwidth_Callback(hObject, eventdata, handles)
% hObject    handle to srcwidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of srcwidth as text
%        str2double(get(hObject,'String')) returns contents of srcwidth as a double


% --- Executes during object creation, after setting all properties.
function srcwidth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to srcwidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in mprof.
function mprof_Callback(hObject, eventdata, handles)
% hObject    handle to mprof (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns mprof contents as cell array
%        contents{get(hObject,'Value')} returns selected item from mprof


% --- Executes during object creation, after setting all properties.
function mprof_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mprof (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mwidth_Callback(hObject, eventdata, handles)
% hObject    handle to mwidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mwidth as text
%        str2double(get(hObject,'String')) returns contents of mwidth as a double


% --- Executes during object creation, after setting all properties.
function mwidth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mwidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in msh_browse.
function msh_browse_Callback(hObject, eventdata, handles)
% hObject    handle to msh_browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[fname fpath] = uigetfile('*.msh','Load mesh file');
set(handles.meshname,'String',fname);


% --- Executes on button press in qm_browse.
function qm_browse_Callback(hObject, eventdata, handles)
% hObject    handle to qm_browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

fname = uigetfile('*.qm','Load QM file');
set(handles.qmname,'String',fname);



function lnampname_Callback(hObject, eventdata, handles)
% hObject    handle to lnampname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lnampname as text
%        str2double(get(hObject,'String')) returns contents of lnampname as a double


% --- Executes during object creation, after setting all properties.
function lnampname_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lnampname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in lnamp_browse.
function lnamp_browse_Callback(hObject, eventdata, handles)
% hObject    handle to lnamp_browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

fname = uigetfile('*.fem','Load data file (log amplitude)');
set(handles.lnampname,'String',fname);



function phasename_Callback(hObject, eventdata, handles)
% hObject    handle to phasename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of phasename as text
%        str2double(get(hObject,'String')) returns contents of phasename as a double


% --- Executes during object creation, after setting all properties.
function phasename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to phasename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in phase_browse.
function phase_browse_Callback(hObject, eventdata, handles)
% hObject    handle to phase_browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

fname = uigetfile('*.fem','Load data file (phase)');
set(handles.phasename,'String',fname);



function freq_Callback(hObject, eventdata, handles)
% hObject    handle to freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of freq as text
%        str2double(get(hObject,'String')) returns contents of freq as a double


% --- Executes during object creation, after setting all properties.
function freq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in runfwd.
function runfwd_Callback(hObject, eventdata, handles)
% hObject    handle to runfwd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

prm = updateprm(handles);
[mdata,pdata] = dotFwd(prm);
showdata(prm,mdata,pdata);


% --- Executes on selection change in linsolver.
function linsolver_Callback(hObject, eventdata, handles)
% hObject    handle to linsolver (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns linsolver contents as cell array
%        contents{get(hObject,'Value')} returns selected item from linsolver


% --- Executes during object creation, after setting all properties.
function linsolver_CreateFcn(hObject, eventdata, handles)
% hObject    handle to linsolver (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function lintol_Callback(hObject, eventdata, handles)
% hObject    handle to lintol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lintol as text
%        str2double(get(hObject,'String')) returns contents of lintol as a double


% --- Executes during object creation, after setting all properties.
function lintol_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lintol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function grid_nx_Callback(hObject, eventdata, handles)
% hObject    handle to grid_nx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of grid_nx as text
%        str2double(get(hObject,'String')) returns contents of grid_nx as a double


% --- Executes during object creation, after setting all properties.
function grid_nx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to grid_nx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function grid_ny_Callback(hObject, eventdata, handles)
% hObject    handle to grid_ny (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of grid_ny as text
%        str2double(get(hObject,'String')) returns contents of grid_ny as a double


% --- Executes during object creation, after setting all properties.
function grid_ny_CreateFcn(hObject, eventdata, handles)
% hObject    handle to grid_ny (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function grid_nz_Callback(hObject, eventdata, handles)
% hObject    handle to grid_nz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of grid_nz as text
%        str2double(get(hObject,'String')) returns contents of grid_nz as a double


% --- Executes during object creation, after setting all properties.
function grid_nz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to grid_nz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes on button press in prm_browse.
function prm_browse_Callback(hObject, eventdata, handles)
% hObject    handle to prm_browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[fname fpath] = uigetfile('*.prm','Load parameter file');
prmfile = [fpath fname];
if strcmpi(fpath, [pwd fpath(length(fpath))]) == 0
    disp('Warning: parameter file not in current working directory.');
    disp(['Changing directory to ' fpath]);
    path(pwd,path); % keep current directory available
    eval(['cd ' ['''' fpath '''']]);
end
prm = toastParam(prmfile);
%prm.prmfile = prmfile;
updategui(handles,prm);
set(handles.figure1,'name',['fwd [' fname ']']);
setappdata(handles.figure1,'prmfile',prmfile);
setappdata(handles.figure1,'prm',prm);


% --- Executes on button press in prm_save.
function prm_save_Callback(hObject, eventdata, handles)
% hObject    handle to prm_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

prmfile = getappdata(handles.figure1,'prmfile');
prm = updateprm(handles);
prm.Write(prmfile);

% --- Executes on button press in prm_saveas.
function prm_saveas_Callback(hObject, eventdata, handles)
% hObject    handle to prm_saveas (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[fname fpath] = uiputfile('*.prm','Save parameter file');
prmfile = [fpath fname];
if strcmpi(fpath, [pwd fpath(length(fpath))]) == 0
    disp('Warning: parameter file not in current working directory.');
    disp(['Changing directory to ' fpath]);
    path(pwd,path); % keep current directory available
    eval(['cd ' ['''' fpath '''']]);
end
prm = updateprm(handles);
prm.Write(prmfile);
set(handles.figure1,'name',['fwd [' fname ']']);
setappdata(handles.figure1,'prmfile',prmfile);


% --- Executes on selection change in prm_mua.
function prm_mua_Callback(hObject, eventdata, handles)
% hObject    handle to prm_mua (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns prm_mua contents as cell array
%        contents{get(hObject,'Value')} returns selected item from prm_mua


% --- Executes during object creation, after setting all properties.
function prm_mua_CreateFcn(hObject, eventdata, handles)
% hObject    handle to prm_mua (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function val_mua_Callback(hObject, eventdata, handles)
% hObject    handle to val_mua (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of val_mua as text
%        str2double(get(hObject,'String')) returns contents of val_mua as a double


% --- Executes during object creation, after setting all properties.
function val_mua_CreateFcn(hObject, eventdata, handles)
% hObject    handle to val_mua (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in prm_mus.
function prm_mus_Callback(hObject, eventdata, handles)
% hObject    handle to prm_mus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns prm_mus contents as cell array
%        contents{get(hObject,'Value')} returns selected item from prm_mus


% --- Executes during object creation, after setting all properties.
function prm_mus_CreateFcn(hObject, eventdata, handles)
% hObject    handle to prm_mus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function val_mus_Callback(hObject, eventdata, handles)
% hObject    handle to val_mus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of val_mus as text
%        str2double(get(hObject,'String')) returns contents of val_mus as a double


% --- Executes during object creation, after setting all properties.
function val_mus_CreateFcn(hObject, eventdata, handles)
% hObject    handle to val_mus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in prm_ref.
function prm_ref_Callback(hObject, eventdata, handles)
% hObject    handle to prm_ref (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns prm_ref contents as cell array
%        contents{get(hObject,'Value')} returns selected item from prm_ref


% --- Executes during object creation, after setting all properties.
function prm_ref_CreateFcn(hObject, eventdata, handles)
% hObject    handle to prm_ref (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function val_ref_Callback(hObject, eventdata, handles)
% hObject    handle to val_ref (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of val_ref as text
%        str2double(get(hObject,'String')) returns contents of val_ref as a double


% --- Executes during object creation, after setting all properties.
function val_ref_CreateFcn(hObject, eventdata, handles)
% hObject    handle to val_ref (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


