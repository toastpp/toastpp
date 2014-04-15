function varargout = recon_gui(varargin)
% RECON_GUI M-file for recon_gui.fig
%      RECON_GUI, by itself, creates a new RECON_GUI or raises the existing
%      singleton*.
%
%      H = RECON_GUI returns the handle to a new RECON_GUI or the handle to
%      the existing singleton*.
%
%      RECON_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RECON_GUI.M with the given input arguments.
%
%      RECON_GUI('Property','Value',...) creates a new RECON_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before recon_gui_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to recon_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help recon_gui

% Last Modified by GUIDE v2.5 12-Jan-2009 12:57:58

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @recon_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @recon_gui_OutputFcn, ...
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


% --- Executes just before recon_gui is made visible.
function recon_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to recon_gui (see VARARGIN)

% Choose default command line output for recon_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

global prm;
prm = toastParam;
if nargin >= 4
    prmfile = varargin{1};
    set(handles.figure1,'name',['recon [' prmfile ']']);
    prm = toastParam(prmfile);
    updategui(handles,prm);
    setappdata(hObject,'prmfile',prmfile);
end
setappdata(hObject,'prm',prm);

% UIWAIT makes recon_gui wait for user response (see UIRESUME)
%uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = recon_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


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
set(handles.listbox1,'String',prmstr(prm,'',''));
% set(handles.meshname,'String',prm.basis.meshfile);
% set(handles.qmname,'String',prm.meas.qmfile);
% set(handles.srctype,'Value',find(strcmpi(get(handles.srctype,'String'),prm.meas.src.type)));
% set(handles.srcprof,'Value',find(strcmpi(get(handles.srcprof,'String'),prm.meas.src.prof)));
% set(handles.srcwidth,'String',prm.meas.src.width);
% set(handles.mprof,'Value',find(strcmpi(get(handles.mprof,'String'),prm.meas.det.prof)));
% set(handles.mwidth,'String',prm.meas.det.width);
% set(handles.linsolver,'Value',find(strcmpi(get(handles.linsolver,'String'),prm.linsolver.method)));
% set(handles.lintol,'String',prm.linsolver.tol);
% set(handles.grid_nx,'String',prm.basis.bdim(1));
% set(handles.grid_ny,'String',prm.basis.bdim(2));
% if length(prm.basis.bdim) >= 3
%     set(handles.grid_nz,'String',prm.basis.bdim(3));
% end
% set(handles.solver,'Value',find(strcmpi({'LM' 'PCG' 'LBFGS'},prm.solver.method)));
% set(handles.solvtol,'String',prm.solver.tol);
% idx = find(strcmpi(get(handles.regul,'String'),prm.regul.method));
% if length(idx) == 0, idx = 1; end
% set(handles.regul,'Value',idx);
% updategui_regul(idx,handles);
% updategui_initprm(prm.initprm.mua, handles.prm_mua, handles.val_mua);
% updategui_initprm(prm.initprm.mus, handles.prm_mus, handles.val_mus);
% updategui_initprm(prm.initprm.ref, handles.prm_ref, handles.val_ref);


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


% Checks for any missing regularisation parameters and fills them
% with default values
function p = regul_initprm(prm)
p = prm;
if strcmpi(p.method,'None') == 0
    if isfield(p,'tau') == 0 || length(p.tau) == 0
        p.tau = 1e-4;
    end
    switch(p.method)
        case 'TV'
            if isfield(p,'tv') == 0
                p.tv = struct ('beta',0.01);
            end
        case 'Huber'
            if isfield(p,'huber') == 0
                p.huber = struct ('eps',0.01);
            end
    end
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


function str=prmstr (prm,prefix,str)
fields = fieldnames(prm);
for i = 1:size(fields,1)
    f = getfield(prm,fields{i});
    if strcmpi(fields{i},'callback'), f = '<internal structure>'; end
    if isstruct(f)
        str = prmstr(f,[prefix fields{i} '.'],str);
    else
        tagstr = sprintf('%25s :', [prefix fields{i}]);
        if ischar(f), valstr = sprintf(' %s', f);
        elseif islogical(f)
            if f == true, valstr = sprintf(' true');
            else, valstr = sprintf(' false'); end
        elseif isreal(f)
            if length(f) > 3
                valstr = sprintf(' [%dx%d real]', size(f,1), size(f,2));
            else
                valstr = sprintf(' %g', f);
            end
        elseif isinteger(f), valstr = sprintf(' %d', f);
        end
        line = [tagstr valstr];
        str{length(str)+1} = line;
    end
end


% --- Executes on button press in runrecon.
function runrecon_Callback(hObject, eventdata, handles)
% hObject    handle to runrecon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
prm = getappdata(handles.figure1,'prm');
prm.transient.callback.iter = @disp_iter;
toastRecon(prm);

function disp_iter(prm,res)
figure(1);
set(gcf,'Name','toastRecon');
bdim = prm.solver.basis.bdim;
if length(bdim) == 2
    subplot(2,2,1), imagesc(reshape(res.bmua,bdim)), axis equal, axis tight, colorbar;
    set(gca,'Units','normalized','OuterPosition',[0 0.5 0.5 0.5]);
    title('initial \mu_a');
    subplot(2,2,2), imagesc(reshape(res.bmus,bdim)), axis equal, axis tight, colorbar;
    set(gca,'Units','normalized','OuterPosition',[0.5 0.5 0.5 0.5]);
    title('initial \mu_s');
    subplot(2,1,2);
    semilogy(res.of);xlabel('iteration');axis tight;title(['Objective:']);
    set(gca,'Units','normalized','OuterPosition',[0 0 1 0.5]);
    drawnow;
else
    muamin=min(res.bmua);
    muamax=max(res.bmua);
    musmin=min(res.bmus);
    musmax=max(res.bmus);
    bmua = reshape(res.bmua,bdim');
    bmus = reshape(res.bmus,bdim');
    for i=1:4
        idx = round(bdim(3)/4*(i-0.5));
        subplot(4,2,i); imagesc(bmua(:,:,idx),[muamin muamax]); axis equal tight
        subplot(4,2,i+4); imagesc(bmus(:,:,idx),[musmin musmax]); axis equal tight
    end
    drawnow
    save recon_gui_res.mat res
end


% --- Executes on button press in prm_browse.
function prm_browse_Callback(hObject, eventdata, handles)
% hObject    handle to prm_browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[fname fpath] = uigetfile('*.prm','Load parameter file');
if ~isequal(fname,0) && ~isequal(fpath,0)
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
    set(handles.figure1,'name',['recon [' fname ']']);
    setappdata(handles.figure1,'prmfile',prmfile);
    setappdata(handles.figure1,'prm',prm);
end


% --- Executes on button press in prm_save.
function prm_save_Callback(hObject, eventdata, handles)
% hObject    handle to prm_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

prmfile = getappdata(handles.figure1,'prmfile');
prm = getappdata(handles.figure1,'prm');
prm.Write(prmfile);


% --- Executes on button press in prm_saveas.
function prm_saveas_Callback(hObject, eventdata, handles)
% hObject    handle to prm_saveas (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[fname fpath] = uiputfile('*.prm','Save parameter file');
if ~isequal(fname,0) && ~isequal(fpath,0)
    prmfile = [fpath fname];
    if strcmpi(fpath, [pwd fpath(length(fpath))]) == 0
        disp('Warning: parameter file not in current working directory.');
        disp(['Changing directory to ' fpath]);
        path(pwd,path); % keep current directory available
        eval(['cd ' ['''' fpath '''']]);
    end
    prm = getappdata(handles.figure1,'prm');
    prm.Write(prmfile);
    set(handles.figure1,'name',['recon [' fname ']']);
    setappdata(handles.figure1,'prmfile',prmfile);
end


% --- Executes on button press in pushbutton13.
function pushbutton13_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global prm
prm = getappdata(handles.figure1,'prm');
assignin('base','prm',prm);
fprintf('Reconstruction parameters exported\nto variable "prm" in workspace.\n');



% =============================================================
% recursively display prm structure fields
function dispprm(prm,prefix)
fields = fieldnames(prm);
for i=1:size(fields,1)
    f = getfield(prm,fields{i});
    if strcmpi(fields{i},'callback'), f = '<internal structure>'; end
    if isstruct(f)
        dispprm(f,[prefix fields{i} '.']);
    else
        fprintf (1, '%30s :', [prefix fields{i}]);
        if     ischar(f),    fprintf (1, ' %s', f);
        elseif islogical(f)
            if f == true, fprintf (1, ' true');
            else, fprintf (1, ' false'); end
        elseif isreal(f)
            if length(f) > 3
                fprintf (1, ' [%dx%d real]', size(f,1), size(f,2));
            else
                fprintf (1, ' %f', f);
            end
        elseif isinteger(f), fprintf (1, ' %d', f);
        end
        fprintf (1, '\n');
    end
end


% --- Executes on button press in pushbutton14.
function pushbutton14_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
prm = getappdata(handles.figure1,'prm');
if isfield(prm,'data')
    tmp = recon_data_gui(prm.data);
else
    tmp = recon_data_gui();
end
if isstruct(tmp)
    prm.data = tmp;
    setappdata(handles.figure1,'prm',prm);
    updategui(handles,prm);
end


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1


% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in pushbutton15.
function pushbutton15_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
prm = getappdata(handles.figure1,'prm');
if isfield(prm,'meas')
    tmp = meas_gui(prm.meas);
else
    tmp = meas_gui();
end
if isstruct(tmp)
    prm.meas = tmp;
    setappdata(handles.figure1,'prm',prm);
    updategui(handles,prm);
end


% --- Executes on button press in pushbutton16.
function pushbutton16_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
prm = getappdata(handles.figure1,'prm');
if isfield(prm,'fwdsolver')
    tmp = fwdsolver_gui(prm.fwdsolver);
else
    tmp = fwdsolver_gui();
end
if isstruct(tmp)
    prm.fwdsolver = tmp;
    setappdata(handles.figure1,'prm',prm);
    updategui(handles,prm);
end


% --- Executes on button press in pushbutton17.
function pushbutton17_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
prm = getappdata(handles.figure1,'prm');
if isfield(prm,'initprm')
    tmp = param_gui(prm.initprm);
else
    tmp = param_gui();
end
if isstruct(tmp)
    prm.initprm = tmp;
    setappdata(handles.figure1,'prm',prm);
    updategui(handles,prm);
end


% --- Executes on button press in pushbutton18.
function pushbutton18_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
prm = getappdata(handles.figure1,'prm');
if isfield(prm,'solver')
    tmp = recon_inv_gui(prm.solver);
else
    tmp = recon_inv_gui();
end
if isstruct(tmp)
    prm.solver = tmp;
    setappdata(handles.figure1,'prm',prm);
    updategui(handles,prm);
end


% --- Executes on button press in help.
function help_Callback(hObject, eventdata, handles)
% hObject    handle to help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
eval('doc recon_gui')


% --- Executes on button press in pushbutton20.
function pushbutton20_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
prm = getappdata(handles.figure1,'prm');
if isfield(prm,'regul')
    tmp = regul_gui(prm.regul);
else
    tmp = regul_gui();
end
if isstruct(tmp)
    prm.regul = tmp;
    setappdata(handles.figure1,'prm',prm);
    updategui(handles,prm);
end


function clbk_iter(arg1,arg2,arg3)
    arg1
