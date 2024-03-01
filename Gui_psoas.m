function varargout = Gui_psoas(varargin)
%% Initialize GUI
% GUI_PSOAS MATLAB code for Gui_psoas.fig
%      GUI_PSOAS, by itself, creates a new GUI_PSOAS or raises the existing
%      singleton*.
%
%      H = GUI_PSOAS returns the handle to a new GUI_PSOAS or the handle to
%      the existing singleton*.
%
%      GUI_PSOAS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_PSOAS.M with the given input arguments.
%
%      GUI_PSOAS('Property','Value',...) creates a new GUI_PSOAS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Gui_psoas_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Gui_psoas_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Gui_psoas

% Last Modified by GUIDE v2.5 19-Feb-2019 09:25:53

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Gui_psoas_OpeningFcn, ...
                   'gui_OutputFcn',  @Gui_psoas_OutputFcn, ...
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

function Gui_psoas_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
guidata(hObject, handles);
function varargout = Gui_psoas_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

%% Load data
function pushbutton1_Callback(hObject, eventdata, handles)  
clc
load_data(handles);



%% Insert_Number_Patient
% function pushbutton5_Callback(hObject, eventdata, handles)
% 
% global num_pat_in num_pat_fin
% num_pat_in=get(handles.edit2, 'string');
% num_pat_in=str2double(num_pat_in);
% num_pat_fin=get(handles.edit3, 'string');
% num_pat_fin=str2double(num_pat_fin);

%% ------
function edit2_Callback(hObject, eventdata, handles)

%% ------
function edit2_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% ------
function edit3_Callback(hObject, eventdata, handles)

%% ------
function edit3_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% %% %% %% %% %% %% %% %% %% %% %%
%% %% %% %% %% FLAGS %% %% %% %% %%
%% %% %% %% %% %% %% %% %% %% %% %%

%% FLAG: Save figures 
function checkbox1_Callback(hObject, eventdata, handles)


%% FLAG: Use active contour for total analysis
function checkbox4_Callback(hObject, eventdata, handles)


%% FLAG: Use edge detection for total analysis
function checkbox5_Callback(hObject, eventdata, handles)




%% %% %% %% %% %% %% %% 
%% Analysis  buttons %% 
%% %% %% %% %% %% %% %% 

%% Centers selection
function pushbutton3_Callback(hObject, eventdata, handles)
global folder_data
% TODO: Cristina commentata parte num_pat_in
% psoas_center_g(folder_data, num_pat_in, num_pat_fin); 
psoas_center_c(folder_data);

%% Edge detection segmentation
function pushbutton2_Callback(hObject, eventdata, handles)
global folder_data num_pat_in num_pat_fin
checksavefig=get(handles.checkbox1,'Value');
 psoas_without_center(folder_data, num_pat_in, num_pat_fin, checksavefig);

%% Active_Contour_3D
function pushbutton6_Callback(hObject, eventdata, handles)
global folder_data num_pat_in num_pat_fin
%  active_contour_3D_dx(folder_data, num_pat_in, num_pat_fin);
%  active_contour_3D_sx(folder_data, num_pat_in, num_pat_fin);
checksavefig=get(handles.checkbox1,'Value');
 active_contour_3D_g(folder_data, num_pat_in, num_pat_fin, checksavefig);
 
%% Mask
function pushbutton9_Callback(hObject, eventdata, handles)
global folder_data num_pat_in num_pat_fin
mask(folder_data, num_pat_in, num_pat_fin);

%% Data analysis
function pushbutton13_Callback(hObject, eventdata, handles)
global folder_data num_pat_in num_pat_fin
data_analysis(folder_data, num_pat_in, num_pat_fin);



%% Pet analysis
function pushbutton15_Callback(hObject, eventdata, handles)
global folder_data num_pat_in num_pat_fin
data_analysisPET(folder_data, num_pat_in, num_pat_fin)

%% Plot Areas in Voxel for each slice
function pushbutton16_Callback(hObject, eventdata, handles)
global folder_data num_pat_in num_pat_fin
plot_areas(folder_data, num_pat_in, num_pat_fin)

%% Save Binary Dicoms
function pushbutton17_Callback(hObject, eventdata, handles)
global folder_data num_pat_in num_pat_fin
% saves 3-colors dicom
savedicom(folder_data, num_pat_in, num_pat_fin );

%% Total Analysis
function pushbutton18_Callback(hObject, eventdata, handles)
global folder_data num_pat_in num_pat_fin
% checks
checkedge = 1.*get(handles.checkbox4, 'Value');
checkactcont = 2.*get(handles.checkbox5, 'Value');
checktot = checkedge + checkactcont;
switch checktot
    case 0
        warning('Neither segmentation method was selected. The default method is edge.')
        checksavefig = get(handles.checkbox1, 'Value');
        % edge
        psoas_without_center(folder_data, num_pat_in, num_pat_fin,checksavefig);
        % active contour
        active_contour_3D_dx(folder_data, num_pat_in, num_pat_fin);
        active_contour_3D_sx(folder_data, num_pat_in, num_pat_fin);
    case 2
        checksavefig = get(handles.checkbox1, 'Value');
        % edge
        psoas_without_center(folder_data, num_pat_in, num_pat_fin,checksavefig);
    case 1
        % active contour
%         active_contour_3D_dx(folder_data, num_pat_in, num_pat_fin);
%         active_contour_3D_sx(folder_data, num_pat_in, num_pat_fin);
        active_contour_3D_g(folder_data, num_pat_in, num_pat_fin);
    case 3
        warning('Both segmentation methods were selected. The default method is edge.')
        checksavefig = get(handles.checkbox1, 'Value');
        % edge
        psoas_without_center(folder_data, num_pat_in, num_pat_fin,checksavefig);
        % active contour
        active_contour_3D_dx(folder_data, num_pat_in, num_pat_fin);
        active_contour_3D_sx(folder_data, num_pat_in, num_pat_fin);
end

% masks
mask(folder_data, num_pat_in, num_pat_fin);
maskPETCT(folder_data, num_pat_in, num_pat_fin);

% data analysis
data_analysis(folder_data, num_pat_in, num_pat_fin);

% data analysis PET
data_analysisPET(folder_data, num_pat_in, num_pat_fin)

% save binary dicoms
savedicom(folder_data, num_pat_in, num_pat_fin );

% signals end of analysis
if num_pat_in == num_pat_fin
fprintf('End CT - PT analysis of patient %d.\n', num_pat_in);
else
fprintf('End CT - PT analysis of patients %d to %d.\n', ...
    num_pat_in,num_pat_fin);
end




%% Mask_pet 
function pushbutton19_Callback(hObject, eventdata, handles)
global folder_data num_pat_in num_pat_fin
maskPETCT(folder_data, num_pat_in, num_pat_fin);


%% Analysis without PT
function pushbutton20_Callback(hObject, eventdata, handles)
global folder_data num_pat_in num_pat_fin
checkedge = 1.*get(handles.checkbox4, 'Value');
checkactcont = 2.*get(handles.checkbox5, 'Value');
checktot = checkedge + checkactcont;
switch checktot
    case 0
        warning('Neither segmentation method was selected. The default method is edge.')
    case 1
        checksavefig = get(handles.checkbox1, 'Value');
        % edge
        psoas_without_center(folder_data, num_pat_in, num_pat_fin,checksavefig);
    case 2
        % active contour
        active_contour_3D_dx(folder_data, num_pat_in, num_pat_fin);
        active_contour_3D_sx(folder_data, num_pat_in, num_pat_fin);
    case 3
        warning('Both segmentation methods were selected. The default method is edge.')
end


% masks
mask(folder_data, num_pat_in, num_pat_fin);

% data analysis
data_analysis(folder_data, num_pat_in, num_pat_fin);

% save binary dicoms
savedicom(folder_data, num_pat_in, num_pat_fin );

% signals end of analysis
if num_pat_in == num_pat_fin
fprintf('End CT analysis of patient %d.\n', num_pat_in);
else
fprintf('End CT analysis of patients %d to %d.\n', ...
    num_pat_in,num_pat_fin);
end






%% CT and PT analysis.
function pushbutton21_Callback(hObject, eventdata, handles)
global folder_data num_pat_in num_pat_fin
% masks
mask(folder_data, num_pat_in, num_pat_fin);
maskPETCT(folder_data, num_pat_in, num_pat_fin);

% data analysis
data_analysis(folder_data, num_pat_in, num_pat_fin);

% data analysis PET
data_analysisPET(folder_data, num_pat_in, num_pat_fin)









%% %% %% %% %% %% %% %% %% %% %% %%
%% %% %%  Close & Clear %% %% %% %% 
%% %% %% %% %% %% %% %% %% %% %% %%
function pushbutton4_Callback(hObject, eventdata, handles )
close all
clear all
clc
