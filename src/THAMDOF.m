function varargout = THAMDOF(varargin)
% Version 2.0
%   - Added video button
%   - Added Collapse SF computation option
%
% Version 2.1
%   - Video with realistic shapes
%
% Version 3.0
%   - Added P-delta effects
%   - Added switch for Strengh limit
%   - Improved computation of Collapse SF
%
% Version 3.1
%   - Added time when completed
%   - Fixed a bug with the length of ug_int
%   - Fixed a bug when using P-delta and Strengh limit in MDOF
%
% Version 3.2
%   - Added the posibility of exporting the interstory velocities also
%
% Version 3.3
%   - Added export displacement time histories option
%
% Version 3.4
%   - Added export shear forces time histories option
%
% Version 3.5
%   - Fixed a bug when zero-padding the records
%
%
% THAMDOF MATLAB code for THAMDOF.fig
%      THAMDOF, by itself, creates a new THAMDOF or raises the existing
%      singleton*.
%
%      H = THAMDOF returns the handle to a new THAMDOF or the handle to
%      the existing singleton*.
%
%      THAMDOF('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in THAMDOF.M with the given input arguments.
%
%      THAMDOF('Property','Value',...) creates a new THAMDOF or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before THAMDOF_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to THAMDOF_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help THAMDOF

% Last Modified by GUIDE v2.5 27-Aug-2015 14:49:54

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @THAMDOF_OpeningFcn, ...
    'gui_OutputFcn',  @THAMDOF_OutputFcn, ...
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


% --- Executes just before THAMDOF is made visible.
function THAMDOF_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to THAMDOF (see VARARGIN)

clc;
global state Building GM Responses

state = 0;

%set(hObject, 'Units', 'Normalized', 'position', [0.05, 0.05, 0.9, 0.9]);
set(handles.figure1,'Units','Pixels','OuterPosition',get(0,'ScreenSize').*[50, 50, 0.9, 0.9])

% Choose default command line output for THAMDOF
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes THAMDOF wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = THAMDOF_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in ExportButton.
function ExportButton_Callback(hObject, eventdata, handles)
% hObject    handle to ExportButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global state Building GM Responses g lengthT forceT

Parameter = get(handles.ExportParameterMenu,'value')-1;

if Parameter == 0
    msgbox('No parameter selected. Try again.','Error');
    return
end

if state < 3
    msgbox('Run analysis first.','Error');
    return
end

[filename pathname] = uiputfile('*.xlsx','Export Results');
if ~ischar(filename) && ~ischar(pathname)   % If no file
    msgbox('No file selected. Try again.','Error');
    return
end

Ngm = length(GM.dt);
Nst = length(Building.W);
DATAmax = zeros(Nst,Ngm);
DATAmin = DATAmax;
DATA = DATAmax;
switch Parameter
    case 1                                              % Displacement
        for i = 1:Ngm
            DATAmax(:,i) = max(Responses.u{i},[],2);
            DATAmin(:,i) = min(Responses.u{i},[],2);
        end
        Title = ['u [' lengthT ']'];
        Sheet = 'u';
    case 2                                              % Interstory Displ.
        for i = 1:Ngm
            DATAmax(:,i) = max(Responses.ID{i},[],2);
            DATAmin(:,i) = min(Responses.ID{i},[],2);
        end
        Title = ['ID [' lengthT ']'];
        Sheet = 'ID';
    case 3                                              % IDR
        for i = 1:Ngm
            DATAmax(:,i) = max(Responses.IDR{i},[],2);
            DATAmin(:,i) = min(Responses.IDR{i},[],2);
        end
        Title = 'IDR [ ]';
        Sheet = 'IDR';
    case 4                                              % Residual ID
        for i = 1:Ngm
            DATA(:,i) = Responses.ID{i}(:,end);
        end
        Title = ['Residual ID [' lengthT ']'];
        Sheet = 'ResID';
    case 5                                              % Residual IDR
        for i = 1:Ngm
            DATA(:,i) = Responses.IDR{i}(:,end);
        end
        Title = 'Residual IDR [ ]';
        Sheet = 'ResIDR';
    case 6                                              % Velocity
        for i = 1:Ngm
            DATAmax(:,i) = max(Responses.v{i},[],2);
            DATAmin(:,i) = min(Responses.v{i},[],2);
        end
        Title = ['v [' lengthT '/sec]'];
        Sheet = 'v';
    case 7                                              % Interstory Velocity
        for i = 1:Ngm
            DATAmax(:,i) = max(Responses.IV{i},[],2);
            DATAmin(:,i) = min(Responses.IV{i},[],2);
        end
        Title = ['Int. vel. [' lengthT '/sec]'];
        Sheet = 'Int. vel';
    case 8                                              % Total Acceleration
        for i = 1:Ngm
            DATAmax(:,i) = max(Responses.a_t{i},[],2);
            DATAmin(:,i) = min(Responses.a_t{i},[],2);
        end
        Title = 'a_t [g]';
        Sheet = 'a_t';
    case 9                                              % Story Restoring Force
        for i = 1:Ngm
            DATAmax(:,i) = max(Responses.fs_st{i},[],2);
            DATAmin(:,i) = min(Responses.fs_st{i},[],2);
        end
        Title = ['fs [' forceT ']'];
        Sheet = 'fs';
    case 10                                              % IDR Time Histories
        st_exp = str2double(inputdlg('Select story: [numeric value]'));
        if isempty(st_exp)
            return
        end
        while isnan(st_exp)
            waitfor(msgbox('Story must be a numeric value. Try again.','Error'))
            st_exp = str2double(inputdlg('Select story: [numeric value]'));
        end
        while Building.Story ~= st_exp
            waitfor(msgbox('Selected story must match with one of the building stories. Try again.','Error'))
            st_exp = str2double(inputdlg('Select story: [numeric value]'));
        end
        
        Sheet = ['TH_IDR_' num2str(st_exp)];
        
        % Maximum number of points
        Ngm_updated = zeros(1,length(GM.dt));
        for igm = 1:Ngm
            Ngm_updated(igm) = length(Responses.IDR{igm}(st_exp,:));
        end
        
        % DATA matrix for writing
        DATA = NaN(max(Ngm_updated),2*length(GM.dt));
        HEADER = cell(1,2*length(GM.dt));
        NPTS = HEADER;
        HEADER2 = HEADER;
        for igm = 1:Ngm
            HEADER(2*igm-1:2*igm) = {GM.Names{igm},GM.Names{igm}};
            NPTS(2*igm-1:2*igm) = {'Npts = ',Ngm_updated(igm)};
            HEADER2(2*igm-1:2*igm) = {'Time [sec]','IDR'};
            DATA(1:Ngm_updated(igm),2*igm-1) = GM.time{igm}';
            DATA(1:Ngm_updated(igm),2*igm) = Responses.IDR{igm}(st_exp,:)';
        end
        
    case 11                                             % Story Restoring Force Time Histories
        st_exp = str2double(inputdlg('Select story: [numeric value]'));
        if isempty(st_exp)
            return
        end
        while isnan(st_exp)
            waitfor(msgbox('Story must be a numeric value. Try again.','Error'))
            st_exp = str2double(inputdlg('Select story: [numeric value]'));
        end
        while Building.Story ~= st_exp
            waitfor(msgbox('Selected story must match with one of the building stories. Try again.','Error'))
            st_exp = str2double(inputdlg('Select story: [numeric value]'));
        end
        
        Sheet = ['TH_fs_' num2str(st_exp)];
        
        % Maximum number of points
        Ngm_updated = zeros(1,length(GM.dt));
        for igm = 1:Ngm
            Ngm_updated(igm) = length(Responses.fs_st{igm}(st_exp,:));
        end
        
        % DATA matrix for writing
        DATA = NaN(max(Ngm_updated),2*length(GM.dt));
        HEADER = cell(1,2*length(GM.dt));
        NPTS = HEADER;
        HEADER2 = HEADER;
        for igm = 1:Ngm
            HEADER(2*igm-1:2*igm) = {GM.Names{igm},GM.Names{igm}};
            NPTS(2*igm-1:2*igm) = {'Npts = ',Ngm_updated(igm)};
            HEADER2(2*igm-1:2*igm) = {'Time [sec]', ['fs [' forceT ']'];};
            DATA(1:Ngm_updated(igm),2*igm-1) = GM.time{igm}';
            DATA(1:Ngm_updated(igm),2*igm) = Responses.fs_st{igm}(st_exp,:)';
        end
        
    case 12     % SF of Collapse
        DATA = Responses.SF_col
        Title = 'Collapse SF';
        Sheet = 'ColSF';
end

h = waitbar(0,'Exporting Data...');
try
    if Parameter ~= 4 && Parameter ~= 5 && Parameter ~= 10 && Parameter ~= 11 && Parameter ~= 12
        xlswrite([pathname filename],{Title},Sheet,'A1:A1')
        xlswrite([pathname filename],GM.Names,Sheet,'B1')
        waitbar(1/4,h,'Exporting Data...')
        xlswrite([pathname filename],DATAmax,Sheet,'B2')
        xlswrite([pathname filename],(1:Nst)',Sheet,'A2')
        waitbar(1/2,h,'Exporting Data...')
        xlswrite([pathname filename],{Title},Sheet,['A' num2str(4+Nst) ':A' num2str(4+Nst)])
        xlswrite([pathname filename],GM.Names,Sheet,['B' num2str(4+Nst)])
        waitbar(3/4,h,'Exporting Data...')
        xlswrite([pathname filename],DATAmin,Sheet,['B' num2str(5+Nst)])
        xlswrite([pathname filename],(1:Nst)',Sheet,['A' num2str(5+Nst)])
        waitbar(1,h,'Exporting Data...')
    elseif Parameter == 10 || Parameter == 11
        xlswrite([pathname filename],HEADER,Sheet,'A1')
        waitbar(1/4,h,'Exporting Data...')
        xlswrite([pathname filename],NPTS,Sheet,'A2')
        waitbar(2/4,h,'Exporting Data...')
        xlswrite([pathname filename],HEADER2,Sheet,'A3')
        waitbar(3/4,h,'Exporting Data...')
        xlswrite([pathname filename],DATA,Sheet,'A4')
        waitbar(1,h,'Exporting Data...')
    elseif Parameter == 12
        xlswrite([pathname filename],{Title},Sheet,'A1:A1')
        waitbar(1/3,h,'Exporting Data...')
        xlswrite([pathname filename],GM.Names,Sheet,'A2')
        waitbar(2/3,h,'Exporting Data...')
        xlswrite([pathname filename],DATA,Sheet,'A3')
        waitbar(1,h,'Exporting Data...')
    else
        xlswrite([pathname filename],{Title},Sheet,'A1:A1')
        waitbar(1/3,h,'Exporting Data...')
        xlswrite([pathname filename],GM.Names,Sheet,'B1')
        waitbar(2/3,h,'Exporting Data...')
        xlswrite([pathname filename],DATA,Sheet,'B2')
        xlswrite([pathname filename],(1:Nst)',Sheet,'A2')
        waitbar(1,h,'Exporting Data...')
    end
    close(h)
    waitfor(msgbox('Done exporting!','Done'));
catch %#ok<CTCH>
    close(h)
    waitfor(msgbox('The file is used by another process. Close it and try again.','Error'));
    return
end

% --- Executes on selection change in ExportParameterMenu.
function ExportParameterMenu_Callback(hObject, eventdata, handles)
% hObject    handle to ExportParameterMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ExportParameterMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ExportParameterMenu



% --- Executes during object creation, after setting all properties.
function ExportParameterMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ExportParameterMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in BuildingButton.
function BuildingButton_Callback(hObject, eventdata, handles)
% hObject    handle to BuildingButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global state Building GM Responses g lengthT forceT

% waitfor(msgbox(['The file must contain 11 columnns: First column: IM '...
%     '(e.g. Pseudo Spectral Accelerations); Second column: Mean Annual'...
%     ' Frequencies of Exceedance.'],'REMEMBER'));

% Open File
[filename pathname] = uigetfile('*.csv','Building Information');
if ~ischar(filename) && ~ischar(pathname)   % If no file
    msgbox('No file selected.','Error');
    set(handles.BuildingText,'String','Error: Must select a file.');
    return
end

LoadedFile = importdata([pathname filename]);
LoadedData = LoadedFile.data;

if size(LoadedData,2) == 12     % Check number of columns
    % Set units
    if get(handles.UnitsMenu,'value') == 1
        lengthT = 'cm';
        forceT = 'tonf';
        g = 980.665;
    else
        lengthT = 'in';
        forceT = 'kip';
        g = 386.09;
    end
    
    % Fill the Building Information
    Building.Story = LoadedData(:,1);
    [Building.Story,ind] = sort(Building.Story);
    Nst = length(Building.Story);
    Building.h = LoadedData(ind,2);
    Building.H = cumsum(Building.h);
    Building.W = LoadedData(ind,3);
    Building.P = LoadedData(ind,4);
    Building.K = LoadedData(ind,5);
    Building.Fy = LoadedData(ind,6);
    Building.as = LoadedData(ind,7);
    Building.ac = LoadedData(ind,8);
    Building.dcdy = LoadedData(ind,9);
    Building.Xi = LoadedData(ind,10);
    Building.do = LoadedData(ind,11);
    Building.Vo = LoadedData(ind,12);
    Building.Fc = Building.Fy.*(1+(Building.dcdy - 1).*Building.as) - Building.P./Building.h.*Building.dcdy.*Building.Fy./Building.K;
    Building.d_col = Building.Fy./Building.K.*(Building.dcdy)+Building.Fc./((-Building.ac+Building.P./(Building.h.*Building.K)).*Building.K);
    
    % Compute modal information
    [Building.T,Building.Gphi] = EigAnalysis(Building.W,Building.K,g);
    
    plotUndeformed(handles);
    
    axes(handles.BuildingAxes);
    
    % Update Mode pop-up menu
    ModeText = cell(Nst+1,1);   StoryText = ModeText;
    ModeText{1} = 'Undeformed Shape';
    StoryText{1} = 'All Stories';
    for i = 1:Nst
        ModeText{i+1} = ['Mode ' num2str(i) '. T = ' num2str(Building.T(i),'%.3f') ' sec ; Xi = ' num2str(Building.Xi(i)*100) '%'];
        StoryText{i+1} = ['Story ' num2str(i)];
    end
    set(handles.ModeMenu,'String',ModeText)
    set(handles.HideCurvesMenu,'String',StoryText)
    
    % Update State
    state = state + 1;
    
    set(handles.BuildingText,'String','Building Loaded!');
    msgbox('Done Loading','Done');
    
    UnitOptions = cellstr(get(handles.UnitsMenu,'String'));
    
    set(handles.UnitsText,'String',UnitOptions{get(handles.UnitsMenu,'Value')})
    set(handles.UnitsMenu,'Visible','off')
    
else
    
    msgbox('Matrix dimensions do not match.','Error');
    set(handles.BuildingText,'String','Error: Try again.');
    
end

% --- Executes on selection change in ModeMenu.
function ModeMenu_Callback(hObject, eventdata, handles)
% hObject    handle to ModeMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ModeMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ModeMenu

Mode = get(hObject,'Value') - 1;

if Mode == 0
    plotUndeformed(handles)
else
    plotMode(handles,Mode)
end


% --- Executes on button press in GMsButton.
function GMsButton_Callback(hObject, eventdata, handles)
% hObject    handle to GMsButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global state Building GM Responses g lengthT forceT

[filename pathname] = uigetfile('*.csv','Ground Motions');
if ~ischar(filename) && ~ischar(pathname)   % If no file
    msgbox('No file selected.','Error');
    set(handles.BuildingText,'String','Error: Must select a file.');
    return
end

LoadedFile = importdata([pathname filename]);
LoadedData = LoadedFile.data;

Ngm = size(LoadedData,2);
GM.Names = LoadedFile.colheaders;
GM.dt = LoadedData(2,:);
GM.SF = LoadedData(3,:);
GM.Acc = cell(1,Ngm);
GM.time = GM.Acc;
GM.Npts = LoadedData(1,:);
for i = 1:size(LoadedData,2)
    Np = GM.Npts(i);
    if Np+3 > size(LoadedData,1)
        GM.Acc{i} = LoadedData(4:end,i);
        GM.Acc{i}(Np) = 0;
    else
        GM.Acc{i} = LoadedData(4:3+Np,i);
    end
    GM.Acc{i}(isnan(GM.Acc{i})) = 0;
    GM.time{i} = 0:GM.dt(i):(GM.dt(i)*(Np-1));
end

% Update State
state = state + 1;

set(handles.GMText,'String','Ground Motions Loaded!');
GMMenuText = ['Select Ground Motion',GM.Names];
set(handles.PlotGMsMenu,'String',GMMenuText)
msgbox(['Done loading ' num2str(Ngm) ' ground motions'],'Done');


% --- Executes on button press in AnalysisButton.
function AnalysisButton_Callback(hObject, eventdata, handles)
% hObject    handle to AnalysisButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global state Building GM Responses g lengthT forceT dSF

if state < 2
    msgbox('First load the Building Information (Input #1) and the Ground Motions (Input #2).'...
        ,'Error')
    return;
end

Ngm = length(GM.dt);
set(handles.ExportParameterMenu,'value',1);
ExportOptions = get(handles.ExportParameterMenu,'String');

tol = 1E-4;
StrengthLimitCheck = get(handles.StrengthLimitCheck,'value');

u = cell(1,Ngm);
v = u; a_t = u; fs_st = u; fs = u; ID = u; IDR = u; IV = u;
h = waitbar(0,['0 of ' num2str(Ngm) ' ground motions completed']...
    ,'Name','Analyzing...');

tic
if get(handles.CollapseSFCheck,'value') == 0    % Just run each ground motion
    
    for igm = 1:Ngm
        [u{igm}, v{igm}, a_t{igm}, fs_st{igm}, fs{igm}, ~, ~, ~, ~, ~, GM.time{igm}, GM.Acc_up{igm}] = ...
            MDOF_Shear_IMK_seismic...
            (Building.h, Building.W, Building.P, Building.K, Building.Xi, GM.Acc{igm}'.*GM.SF(igm), GM.dt(igm), ...
            Building.do, Building.Vo, Building.Fy, Building.as, Building.dcdy,...
            Building.ac, tol, g, GM.Names{igm},StrengthLimitCheck);
        
        Uaux = [zeros(1,size(u{igm},2));u{igm}];
        ID{igm} = diff(Uaux,[],1);
        IDR{igm} = (1./Building.h)*ones(1,size(ID{igm},2)).*ID{igm};
        Vaux = [zeros(1,size(v{igm},2));v{igm}];
        IV{igm} = diff(Vaux,[],1);
        waitbar(igm/Ngm,h,[num2str(igm) ' of ' num2str(Ngm) ' ground motions completed'])
    end
    close(h);
    
    Responses.u = u;
    Responses.v = v;
    Responses.a_t = a_t;
    Responses.fs_st = fs_st;
    Responses.fs = fs;
    Responses.ID = ID;
    Responses.IDR = IDR;
    Responses.IV = IV;
    
    % Update Menu of Ground Motions, with their SF
    GMMenuText = ['Select Ground Motion',GM.Names];
    for i = 2:length(GMMenuText)
        GMMenuText{i} = [GMMenuText{i} ' (SF = ' num2str(GM.SF(i-1)) ')'];
    end
    
    set(handles.PlotGMsMenu,'String',GMMenuText)
    
    if length(ExportOptions) == 13
        set(handles.ExportParameterMenu,'String',ExportOptions(1:12));
    end
    
else            % Run each ground motion up to collapse
    pos = get(h,'position');
    set(h,'Position',[pos(1) pos(2) pos(3) pos(4)*1.3])
    
    SF = dSF*ones(1,Ngm);
    
    for igm = 1:Ngm     % for each ground motion
        cont = true;
        dSF_use = dSF;
        while cont
            waitbar((igm-1)/Ngm,h,{['Analizing "' GM.Names{igm} '" - SF = ' num2str(SF(igm))]...
                ,[num2str(igm-1) ' of ' num2str(Ngm) ' ground motions completed']})
            ugSF = GM.Acc{igm}'.*SF(igm);
            [u{igm}, v{igm}, a_t{igm}, fs_st{igm}, fs{igm}, ~, ~, ~, ~, ~, GM.time{igm}, GM.Acc_up{igm}] = ...
                MDOF_Shear_IMK_seismic...
                (Building.h, Building.W, Building.P, Building.K, Building.Xi, ugSF, GM.dt(igm), ...
                Building.do, Building.Vo, Building.Fy, Building.as, Building.dcdy,...
                Building.ac, tol, g, GM.Names{igm},StrengthLimitCheck);
            
            Uaux = [zeros(1,size(u{igm},2));u{igm}];
            ID{igm} = diff(Uaux,[],1);
            
            MAX_ID = max(abs(ID{igm}),[],2);
            if MAX_ID <= Building.d_col     % No collapse yet
                SF(igm) = SF(igm) + dSF_use;
            else                            % Collapsed
                if dSF_use/SF(igm) <= 0.01  % Convergence criteria: 1% of difference
                    cont = false;
                else
                    dSF_use = dSF_use/2;
                    SF(igm) = SF(igm) - dSF_use;
                end
            end
            
        end
        IDR{igm} = (1./Building.h)*ones(1,size(ID{igm},2)).*ID{igm};
        Vaux = [zeros(1,size(v{igm},2));v{igm}];
        IV{igm} = diff(Vaux,[],1);
        waitbar(igm/Ngm,h,[num2str(igm) ' of ' num2str(Ngm) ' ground motions completed'])
    end
    close(h);
    
    
    Responses.u = u;
    Responses.v = v;
    Responses.a_t = a_t;
    Responses.fs_st = fs_st;
    Responses.fs = fs;
    Responses.ID = ID;
    Responses.IDR = IDR;
    Responses.IV = IV;
    Responses.SF_col = SF;
    
    % Update Menu of Ground Motions, with their SF
    GMMenuText = ['Select Ground Motion',GM.Names];
    for i = 2:length(GMMenuText)
        GMMenuText{i} = [GMMenuText{i} ' (SF_col = ' num2str(SF(i-1)) ')'];
    end
    
    set(handles.PlotGMsMenu,'String',GMMenuText)
    
    if length(ExportOptions) == 12
        ExportOptions{end+1} = 'Collapse SF';
        set(handles.ExportParameterMenu,'String',ExportOptions);
    end
    
end

AnalysisTime=toc;
% Update State
state = 3;

set(handles.RunText,'String','Analysis Complete!');

msgbox(['Analysis completed in ' num2str(AnalysisTime,'%.1f') ' seconds!'],'Done');

% --- Executes on button press in ExitButton.
function ExitButton_Callback(hObject, eventdata, handles)
% hObject    handle to ExitButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close all;

% --- Executes on selection change in PlotGMsMenu.
function PlotGMsMenu_Callback(hObject, eventdata, handles)
% hObject    handle to PlotGMsMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns PlotGMsMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from PlotGMsMenu

plotResponses(handles)


% --- Executes during object creation, after setting all properties.
function PlotGMsMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PlotGMsMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in PlotParameterYMenu.
function PlotParameterYMenu_Callback(hObject, eventdata, handles)
% hObject    handle to PlotParameterYMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns PlotParameterYMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from PlotParameterYMenu

plotResponses(handles)

if get(hObject,'Value') == 1
    set(handles.YUnits,'string','');
end


% --- Executes during object creation, after setting all properties.
function PlotParameterYMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PlotParameterYMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in PlotParameterXMenu.
function PlotParameterXMenu_Callback(hObject, eventdata, handles)
% hObject    handle to PlotParameterXMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns PlotParameterXMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from PlotParameterXMenu

plotResponses(handles)

% --- Executes during object creation, after setting all properties.
function PlotParameterXMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PlotParameterXMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function PlotStoryMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PlotStoryMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in UnitsMenu.
function UnitsMenu_Callback(hObject, eventdata, handles)
% hObject    handle to UnitsMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns UnitsMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from UnitsMenu


% --- Executes during object creation, after setting all properties.
function UnitsMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to UnitsMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes during object creation, after setting all properties.
function ModeMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ModeMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in HideCurvesMenu.
function HideCurvesMenu_Callback(hObject, eventdata, handles)
% hObject    handle to HideCurvesMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns HideCurvesMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from HideCurvesMenu
HideCurves(handles)

% --- Executes during object creation, after setting all properties.
function HideCurvesMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to HideCurvesMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in .
function CollapseSFCheck_Callback(hObject, eventdata, handles)
% hObject    handle to CollapseSFCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CollapseSFCheck

global dSF

if get(hObject,'Value') == 1
    dSF = str2double(inputdlg('Increment of Scale Factor (dSF): [numeric value]'));
    while isnan(dSF)
        waitfor(msgbox('Increment must be a numeric value. Try again.','Error'))
        dSF = str2double(inputdlg('Increment of Scale Factor (dSF): [numeric value]'));
    end
    while dSF < 0.02
        waitfor(msgbox('Choose a value greater than 0.02','Error'))
        dSF = str2double(inputdlg('Increment of Scale Factor (dSF): [numeric value]'));
    end
    set(handles.dSFText,'string',['dSF = ' num2str(dSF)])
else
    set(handles.dSFText,'string','')
end

% --- Executes on button press in StrengthLimitCheck.
function StrengthLimitCheck_Callback(hObject, eventdata, handles)
% hObject    handle to StrengthLimitCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of StrengthLimitCheck

if get(hObject,'Value') == 0
    waitfor(msgbox({'It is recommended to include the Strength Limit option.',...
        'Only uncheck this checkbox if you are comparing against IIIDAP.'},'Warning'))
end

% --- Executes on button press in VideoButton.
function VideoButton_Callback(hObject, eventdata, handles)
% hObject    handle to VideoButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global state Building GM Responses g lengthT forceT leg

indGM = get(handles.PlotGMsMenu,'value')-1;

if state < 3
    waitfor(msgbox('Run the analysis first.','Error'))
    return
end

if indGM == 0
    waitfor(msgbox('Choose a ground motion.','Error'))
    return
end

indC = ones(length(Building.W),1)*length(GM.time{indGM});
Collapse = false;
% See if the GM produced collapse
for i = 1:length(indC)
    indAux = find(abs(Responses.ID{indGM}(i,:)) >= Building.d_col(i),1);
    if ~isempty(indAux)
        indC(i) = indAux;
        Collapse = true;
    end
end
[indC,StCol] = min(indC);

Nextra = 250;
if indC < length(GM.time{indGM})-Nextra
    indC = indC + Nextra;
else
    Nextra = length(GM.time{indGM}) - indC;
    indC = length(GM.time{indGM});
end

AR = inputdlg('Enter Aspect Ratio for horizontal axes.','Aspect Ratio',1,{'20'});
AR = str2double(AR);

if ~isempty(AR)
    Nst = length(Building.W);
    Xf = zeros(Nst+1,indC);
    Yf = zeros(Nst+1,indC);
    for i = 2:Nst+1
        Xf(i,:) = Responses.u{indGM}(i-1,1:indC);
        Yf(i,:) = Yf(i-1,:) + sqrt(Building.h(i-1)^2 - (Xf(i-1,:)-Responses.u{indGM}(i-1,1:indC)).^2);
    end
    
    ID = Responses.ID{indGM}(:,1:indC)';
    ug = GM.Acc_up{indGM}';
    time = GM.time{indGM}';
    
    l = cell(Nst,1);
    for i = 1:Nst
        l{i} = ['Story ' num2str(i)];
    end
    
    figure('Position',[420   100   800   850])
    a1 = subplot('Position',[0.25 0.46 0.6 0.49]);
    hold on
    for i = 1:Nst
        [Xp,Yp] = PlotIDShape(Xf(i:i+1,1),Yf(i:i+1,1));
        h1_shape(i) = plot(Xp,Yp,'b-','linewidth',3);
        h1_points(i) = plot(Xf(i:i+1,1),Yf(i:i+1,1),'bo','linewidth',3);
    end
    Ylimits = get(a1,'ylim');
    grid on
    xlabel(['Displacement [' lengthT ']'])
    ylabel(['Height [' lengthT ']'])
    ylim(Ylimits)
    title('t = 0.000 sec')
    daspect(a1,[1 AR 1])
    xlim([-1 1]*max(max(abs(Responses.u{indGM}(:,1:indC-Nextra))))*1.5)
    
    
    a2 = subplot('Position',[0.1 0.25 0.85 0.15]);
    h2 = plot(time(1),ID(1,:));
    hold on
    ylabel(['ID [' lengthT ']'])
    grid on
    legend(l,'Location','NorthEast','autoupdate','off')
    xlim([0 1.3*max(time)])
    ylim([min(min((ID(1:(indC-Nextra),:)))) max(max((ID(1:(indC-Nextra),:))))]*1.3)
    
    subplot('Position',[0.1 0.06 0.85 0.15])
    plot(time,ug,'k')
    hold on
    h3 = plot(time(1),ug(1),'r');
    ylabel('Ground Acc. [g]')
    xlabel('Time [sec]')
    grid on
    xlim([0 1.3*max(time)])
    
    for i = 2:size(Xf,2)
        pause((time(i)-time(i-1))/2)
        title(a1,['t = ' num2str(time(i),'%.3f') ' sec'])
        
        for st = 1:Nst
            maxID = max(abs(ID(1:i,st)));
            if maxID <= Building.Fy(st)/Building.K(st)
                ColorState = 'b';
                Line = '-';
            elseif maxID <= Building.Fy(st)/Building.K(st)*Building.dcdy
                ColorState = [1 0.7 0.5];
                Line = '-';
            elseif maxID <= Building.d_col(st)
                ColorState = 'r';
                Line = '-';
            else
                ColorState = 'r';
                Line = '--';
            end
            
            [Xp,Yp] = PlotIDShape(Xf(st:st+1,i),Yf(st:st+1,i));
            
            set(h1_shape(st),'Xdata',Xp);
            set(h1_shape(st),'Ydata',real(Yp));
            set(h1_shape(st),'Color',ColorState);
            set(h1_shape(st),'LineStyle',Line);
            set(h1_points(st),'Xdata',Xf(st:st+1,i));
            set(h1_points(st),'Ydata',real(Yf(st:st+1,i)));
            set(h1_points(st),'Color',ColorState);
            
        end
        
        set(h2,'Xdata',time(1:i));
        for st = 1:length(Building.W)
            set(h2(st),'Ydata',ID(1:i,st));
        end
        set(h3,'Xdata',time(1:i));
        set(h3,'Ydata',ug(1:i));
        
    end
    
    if Collapse
        plot(a2,time(indC-Nextra),ID(indC-Nextra,StCol),'rx','linewidth',2 ...
            ,'MarkerSize',8);
        axes(a2);
        text(time(indC-Nextra),ID(indC-Nextra,StCol),'   Collapse','color','r'...
            ,'FontWeight','bold');
        
        for st = 1:Nst
            i = indC-Nextra;
            maxID = max(abs(ID(1:i,st)));
            if maxID <= Building.Fy(st)/Building.K(st)
                ColorState = 'b';
                Line = '-';
            elseif maxID <= Building.Fy(st)/Building.K(st)*Building.dcdy
                ColorState = [1 0.7 0.5];
                Line = '-';
            elseif maxID <= Building.d_col(st)
                ColorState = 'r';
                Line = '-';
            else
                ColorState = 'r';
                Line = '--';
            end
            
            [Xp,Yp] = PlotIDShape(Xf(st:st+1,i),Yf(st:st+1,i));
            plot(a1,Xp,Yp,'Color',ColorState,'LineStyle',Line,'linewidth',1.5);
            plot(a1,Xf(st:st+1,i),Yf(st:st+1,i),'o','Color',ColorState,'linewidth',1.5)
        end
    end
end


function [Xp,Yp] = PlotIDShape(Xf,Yf)

Np = 20;
Amp = (Xf(1)-Xf(2))/2;
T = (Yf(2)-Yf(1))*2;
Xp = Xf(1)+Amp*cos((2*pi)/(2*Np)*(0:Np))-Amp;
Yp = Yf(1)+T/(2*Np)*(0:Np);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [T,Gphi] = EigAnalysis(wi,k,g)
Nst = length(wi);

% M and K matrices
M = diag(wi)/g;         % Mass Matrix
k_aux = k(2:end);   k_aux(end+1) = 0;
K = diag(k+k_aux) - diag(k(2:end),1) - diag(k(2:end),-1);

[phi,w2] = eig(K,M);
w = sqrt(diag(w2));     % Undamped frequencies
[w,index] = sort(w);
T = 2*pi./w;            % Undamped periods

% Sort vectors (modal shapes) and normalize them at roof: phi_roof = 1.0
sphi = phi;
for i = 1:length(wi)
    sphi(:,i) = phi(:,index(i))/ phi(end,index(i));
end
phi = sphi;             % Normalized modal shapes

r = ones(Nst,1);
Gamma = phi'*M*r./diag(phi'*M*phi);

Gphi = phi;
for i = 1:Nst
    Gphi(:,i) = Gamma(i)*phi(:,i);
end


function [] = plotUndeformed(handles)
global Building g lengthT forceT

Nst = length(Building.Story);

plot(handles.BuildingAxes,[0;Building.do],[0;Building.H],'-o','linewidth',2)
axes(handles.BuildingAxes);
grid on
xlabel(['do [' lengthT ']'])
ylabel(['Height [' lengthT ']'])
ylim([0 1.2*Building.H(end)])

doAux = [0;Building.do];
HAux = [0;Building.H];
for i = 1:Nst
    Xtext = (doAux(i) + doAux(i+1))/2;
    Ytext = (HAux(i) + HAux(i+1))/2;
    text(Xtext,Ytext,['   K_' num2str(i) ' = ' num2str(Building.K(i)) ' ' forceT '/' lengthT])
    Xtext = Building.do(i);
    Ytext = Building.H(i);
    text(Xtext,Ytext,['   W_' num2str(i) ' = ' num2str(Building.W(i)) ' ' forceT])
end


function [] = plotMode(handles,Mode)
global Building g lengthT forceT

plot(handles.BuildingAxes,[0;Building.Gphi(:,Mode)],[0;Building.H],'-o','linewidth',2)
axes(handles.BuildingAxes);
grid on
xlabel(['\Gamma_' num2str(Mode) ' \phi_' num2str(Mode)])
ylabel(['Height [' lengthT ']'])
ylim([0 1.2*Building.H(end)])
xlim(1.2*[-max(max(abs(Building.Gphi))) max(max(abs(Building.Gphi)))]);



function [] = plotResponses(handles)
global state Building GM Responses g lengthT forceT leg

ParameterX = get(handles.PlotParameterXMenu,'value');
ParameterY = get(handles.PlotParameterYMenu,'value');
indGM = get(handles.PlotGMsMenu,'value')-1;

if indGM == 0 || ParameterY == 1
    cla(handles.ResultsAxes)
    legend off
    return
end

if ParameterY > 2 && state < 3
    msgbox('First analize the building.','Error')
    set(handles.PlotParameterYMenu,'value',1)
    set(handles.YUnits,'String','')
    return
end

legTot = cell(1,length(Building.W));
for st = 1:length(Building.W)
    legTot{st} = ['Story ' num2str(st)];
end
leg = legTot;


indC = ones(length(Building.W),1)*length(GM.time{indGM});
if state == 3   % See if the GM produced collapse
    for i = 1:length(indC)
        indAux = find(abs(Responses.ID{indGM}(i,:)) >= Building.d_col(i),1);
        if ~isempty(indAux)
            indC(i) = indAux;
        end
    end
end

switch ParameterY
    case 2  % Ground Motion
        if state < 3
            PlotY = GM.Acc{indGM}';
            indC = length(PlotY);
        else
            PlotY = GM.Acc_up{indGM}';
            indC = length(PlotY);
        end
        Ytext = '[g]';
        leg = 'Ground Acceleration';
        
    case 3  % Displacement
        PlotY = Responses.u{indGM}(:,1:min(indC))';
        Ytext = ['[' lengthT ']'];
        
    case 4  % Interstory disp
        PlotY = Responses.ID{indGM}(:,1:min(indC))';
        Ytext = ['[' lengthT ']'];
        
    case 5  % IDR
        PlotY = Responses.IDR{indGM}(:,1:min(indC))';
        Ytext = '[ ]';
        
    case 6  % Velocity
        PlotY = Responses.v{indGM}(:,1:min(indC))';
        Ytext = ['[' lengthT '/sec]'];
        
    case 7  % Total Acceleration
        PlotY = Responses.a_t{indGM}(:,1:min(indC))';
        Ytext = '[g]';
        
    case 8  % Restoring Force
        PlotY = Responses.fs_st{indGM}(:,1:min(indC))';
        Ytext = ['[' forceT ']'];
        
    case 9  % Base Shear
        PlotY = Responses.fs_st{indGM}(1,1:min(indC))';
        Ytext = ['[' forceT ']'];
        leg = 'Base Shear';
        
end

switch ParameterX
    case 1  % Time
        PlotX = GM.time{indGM}(:,1:min(indC))';
        Xtext = '[sec]';
        
    case 2  % Interstory Displacement
        PlotX = Responses.ID{indGM}(:,1:min(indC))';
        Xtext = ['[' lengthT ']'];
        leg = legTot;
        
    case 3  % IDR
        PlotX = Responses.IDR{indGM}(:,1:min(indC))';
        Xtext = '[ ]';
        leg = legTot;
        
end

if ParameterX == 1
    xlimit = [0 max(GM.time{indGM})*1.2];
else
    xlimit = [min(min(PlotX)) max(max(PlotX))]*1.2;
end
ylimit = [min(min(PlotY)) max(max(PlotY))]*1.2;

plot(handles.ResultsAxes,PlotX,PlotY,'linewidth',1);
axes(handles.ResultsAxes);
grid on
legend(leg,'Location','best')
set(handles.XUnits,'String',Xtext)
set(handles.YUnits,'String',Ytext);
xlim(xlimit)
ylim(ylimit)

HideCurves(handles)


function [] = HideCurves(handles)
global leg Building

Nst = length(Building.W);

if isempty(get(handles.ResultsAxes,'children'))
    return
end

Story = get(handles.HideCurvesMenu,'value')-1;
Curves = get(handles.ResultsAxes,'children');

if Story == 0 || length(Curves) == 1
    set(Curves,'visible','on')
    legend(leg,'Location','Best')
else
    set(Curves,'visible','off')
    set(Curves(Nst+1-Story),'visible','on')
    legend(Curves(Nst+1-Story),['Story ' num2str(Story)],'Location','Best')
end
