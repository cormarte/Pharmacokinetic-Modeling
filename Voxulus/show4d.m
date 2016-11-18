function varargout = show4d(varargin)
% SHOW4D M-file for show4d.fig
%      SHOW4D, by itself, creates a new SHOW4D or raises the existing
%      singleton*.
%
%      H = SHOW4D returns the handle to a new SHOW4D or the handle to
%      the existing singleton*.
%
%      SHOW4D('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SHOW4D.M with the given input arguments.
%
%      SHOW4D('Property','Value',...) creates a new SHOW4D or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before show4d_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to show4d_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help show4d

% Last Modified by GUIDE v2.5 01-Mar-2005 13:15:42

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @show4d_OpeningFcn, ...
                   'gui_OutputFcn',  @show4d_OutputFcn, ...
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

function updateAxes (handles)
  set (gcf, 'Name', sprintf ('show4d: t=%.0f z=%.0f y=%.0f x=%.0f', handles.t, handles.z, handles.y, handles.x));
  nx = size(handles.data,3);
  ny = size(handles.data,2);
  nz = size(handles.data,1);
  nt = size(handles.data,4);
  set (handles.xSlider, 'Min', 1);
  set (handles.xSlider, 'Max', size (handles.data, 3));
  set (handles.xSlider, 'SliderStep', [ 1.0/(nx-1) 1 ]);
  set (handles.xSlider, 'Value', handles.x);
  set (handles.ySlider, 'Min', 1);
  set (handles.ySlider, 'Max', size (handles.data, 2));
  set (handles.ySlider, 'SliderStep', [ 1.0/(ny-1) 1 ]);
  set (handles.ySlider, 'Value', handles.y);
  set (handles.zSlider, 'Min', 1);
  set (handles.zSlider, 'Max', size (handles.data, 1));
  set (handles.zSlider, 'SliderStep', [ 1.0/(nz-1) 1 ]);
  set (handles.zSlider, 'Value', handles.z);
  set (handles.tSlider, 'Min', 1);
  set (handles.tSlider, 'Max', size (handles.data, 4));
  set (handles.tSlider, 'SliderStep', [ 1.0/(nt-1) 1 ]);
  set (handles.tSlider, 'Value', handles.t);
  if get (handles.viewModeRadiobutton, 'Value')
    viewMin = min (handles.data(:));
    set (handles.viewMinSlider, 'Min', viewMin);
    set (handles.viewMaxSlider, 'Min', viewMin);
    set (handles.viewMinSlider, 'Value', viewMin);
    viewMax = max (handles.data(:));
    set (handles.viewMinSlider, 'Max', viewMax);
    set (handles.viewMaxSlider, 'Max', viewMax);
    set (handles.viewMaxSlider, 'Value', viewMax);
  else
    viewMin = get (handles.viewMinSlider, 'Value');
    viewMax = get (handles.viewMaxSlider, 'Value');
  end
  d = squeeze (handles.data(handles.z,:,:,handles.t));
  axes (handles.colAxes);
  img = imagesc (d, [ viewMin viewMax ]);
  set (gca, 'Visible', 'off');
  set (img, 'Visible', 'off');
  colorbar ('WestOutside');
  axes (handles.xyAxes);
  v = logical (zeros (size (d)));
  if get (handles.voiRadiobutton, 'Value') && ~isempty (handles.voi)
    v = squeeze (handles.voi(handles.z,:,:));
  else
    v(handles.y,:) = 1;
    v(:,handles.x) = 1;
  end
  show_image (d, [ viewMin viewMax ], v);
  axis xy;
  set (get (handles.xyAxes, 'Children'), 'ButtonDownFcn', {@xyAxes_ButtonDownFcn,handles});
  axes (handles.yzAxes);
  d = squeeze (handles.data(:,:,handles.x,handles.t));
  v = logical (zeros (size (d)));
  if get (handles.voiRadiobutton, 'Value') && ~isempty (handles.voi)
    v = squeeze (handles.voi(:,:,handles.x));
  else
    v(handles.z,:) = 1;
    v(:,handles.y) = 1;
  end
  show_image (d', [ viewMin viewMax ], v');
  axis xy;
  set (get (handles.yzAxes, 'Children'), 'ButtonDownFcn', {@yzAxes_ButtonDownFcn,handles});
  axes (handles.xzAxes);
  d = squeeze (handles.data(:,handles.y,:,handles.t));
  v = logical (zeros (size (d)));
  if get (handles.voiRadiobutton, 'Value') && ~isempty (handles.voi)
    v = squeeze (handles.voi(:,handles.y,:));
  else
    v(handles.z,:) = 1;
    v(:,handles.x) = 1;
  end
  show_image (d, [ viewMin viewMax ], v);
  set (get (handles.xzAxes, 'Children'), 'ButtonDownFcn', {@xzAxes_ButtonDownFcn,handles});
  axes (handles.tacAxes);
  d = squeeze(handles.data(handles.z,handles.y,handles.x,:));
  if get (handles.fitRadiobutton, 'Value')
    t = get (handles.voxDataPopupmenu, 'String');
    dataName = t{get (handles.voxDataPopupmenu, 'Value')};
    t = get (handles.voxErrorPopupmenu, 'String');
    errorName = t{get (handles.voxErrorPopupmenu, 'Value')};
    tac = voxulus ('get', 'tac', dataName, errorName, [ handles.x handles.y handles.z ]);
    c = voxulus ('create', 'tac', dataName, [ handles.x handles.y handles.z ]);
    [ nt nc ] = size (c);
    t = tac(:,1);
    d = [ d, tac(:,2) c(:,3:nc) ];
    plot (t, d);
    axis ([ t(1) t(end) viewMin viewMax ]);
    leg{1} = 'data';
    leg{2} = 'fit';
    for i=3:nc-1
      leg{i} = sprintf ('c%d(t)', i-2);
    end
    leg{nc} = 'cp(t)';
    legend (leg);
  else
    t = 1:size(handles.data,4);
    plot (t, d);
    axis ([ 1 size(handles.data,4) viewMin viewMax ]);
  end
  title (sprintf ('TAC x=%.0f y=%.0f z=%.0f', handles.x, handles.y, handles.z));


% --- Executes just before show4d is made visible.
function show4d_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to show4d (see VARARGIN)

% Choose default command line output for show4d
handles.output = hObject;
handles.data = varargin{1};
if length (varargin) > 1
  handles.voi = varargin{2};
else
  handles.voi = [];
end
handles.dataMin = min (handles.data(:));
handles.dataMax = max (handles.data(:));
set (handles.viewMinSlider, 'Min', handles.dataMin);
set (handles.viewMinSlider, 'Max', handles.dataMax);
set (handles.viewMinSlider, 'Value', handles.dataMin);
set (handles.viewMaxSlider, 'Max', handles.dataMax);
set (handles.viewMaxSlider, 'Min', handles.dataMin);
set (handles.viewMaxSlider, 'Value', handles.dataMax);  
nx = size(handles.data,3);
ny = size(handles.data,2);
nz = size(handles.data,1);
nt = size(handles.data,4);
handles.t = 1;
handles.x = round (nx / 2);
handles.y = round (ny / 2);
handles.z = round (nz / 2);
handles.voxData = '';
handles.voxError = '';
t = voxulus ('get', 'data');
t{length(t)+1} = 'none';
set (handles.voxDataPopupmenu, 'String', t);
set (handles.voxDataPopupmenu, 'Value', length(t));
set (handles.voxErrorPopupmenu, 'String', t);
set (handles.voxErrorPopupmenu, 'Value', length(t));
updateAxes (handles);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes show4d wait for user response (see UIRESUME)
%uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = show4d_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
%varargout{1} = handles.data;


% --- Executes on slider movement.
function xSlider_Callback(hObject, eventdata, handles)
% hObject    handle to xSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
  handles.x = round (get (hObject,'Value'));
  updateAxes (handles);
  guidata(hObject, handles);


% --- Executes on slider movement.
function ySlider_Callback(hObject, eventdata, handles)
% hObject    handle to ySlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
  handles.y = round (get (hObject,'Value'));
  updateAxes (handles);
  guidata(hObject, handles);


% --- Executes on slider movement.
function zSlider_Callback(hObject, eventdata, handles)
% hObject    handle to zSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
  handles.z = round (get (hObject,'Value'));
  updateAxes (handles);
  guidata(hObject, handles);


% --- Executes on selection change in colormapPopupmenu.
function colormapPopupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to colormapPopupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns colormapPopupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from colormapPopupmenu
  t = get (hObject, 'String');
  map = t{get (hObject, 'Value')};
  if strcmp (map, 'gray')
    colormap (gray);
  elseif strcmp (map, 'hot')
    colormap (hot);
  elseif strcmp (map, 'hsv')
    colormap (hsv);
  elseif strcmp (map, 'cool')
    colormap (cool);
  elseif strcmp (map, 'bone')
    colormap (bone);
  elseif strcmp (map, 'copper')
    colormap (copper);
  elseif strcmp (map, 'jet')
    colormap (jet);
  elseif strcmp (map, 'pink')
    colormap (pink);
  elseif strcmp (map, 'autumn')
    colormap (autumn);
  elseif strcmp (map, 'spring')
    colormap (spring);
  elseif strcmp (map, 'summer')
    colormap (summer);
  elseif strcmp (map, 'winter')
    colormap (winter);
  end
  updateAxes (handles);
  guidata (gcf, handles);


% --- Executes on mouse press over axes background.
function xyAxes_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to xyAxes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  p = get (gca, 'CurrentPoint');
  x = round (p(1,1));
  y = round (p(1,2));
  handles.x = x;
  handles.y = y;
  updateAxes (handles);
  guidata(gcf, handles);


% --- Executes on mouse press over axes background.
function yzAxes_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to yzAxes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  p = get (gca, 'CurrentPoint');
  z = round (p(1,1));
  y = round (p(1,2));
  handles.z = z;
  handles.y = y;
  updateAxes (handles);
  guidata (gcf, handles);


% --- Executes on mouse press over axes background.
function xzAxes_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to yzAxes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  p = get (gca, 'CurrentPoint');
  x = round (p(1,1));
  z = round (p(1,2));
  handles.x = x;
  handles.z = z;
  updateAxes (handles);
  guidata (gcf, handles);


% --- Executes on slider movement.
function viewMaxSlider_Callback(hObject, eventdata, handles)
% hObject    handle to viewMaxSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
  set (handles.viewMinSlider, 'Max', get (hObject,'Value'));
  updateAxes (handles);
  guidata(hObject, handles);


% --- Executes on slider movement.
function viewMinSlider_Callback(hObject, eventdata, handles)
% hObject    handle to viewMinSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
  set (handles.viewMaxSlider, 'Min', get (hObject,'Value'));
  updateAxes (handles);
  guidata(hObject, handles);


% --- Executes on button press in viewModeRadiobutton.
function viewModeRadiobutton_Callback(hObject, eventdata, handles)
% hObject    handle to viewModeRadiobutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of viewModeRadiobutton
  if get (hObject, 'Value')
    set (handles.viewMinSlider, 'Enable', 'inactive');
    set (handles.viewMaxSlider, 'Enable', 'inactive');
  else
    set (handles.viewMinSlider, 'Enable', 'on');
    set (handles.viewMaxSlider, 'Enable', 'on');
  end
  updateAxes (handles);
  guidata(hObject, handles);


% --- Executes on slider movement.
function tSlider_Callback(hObject, eventdata, handles)
% hObject    handle to tSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
  handles.t = round (get (hObject,'Value'));
  updateAxes (handles);
  guidata(hObject, handles);


% --- Executes on button press in undoPushbutton.
function undoPushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to undoPushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  load show4d_data data;
  handles.data = data;
    handles.dataMin = min (handles.data(:));
    handles.dataMax = max (handles.data(:));
    set (handles.viewMinSlider, 'Min', handles.dataMin);
    set (handles.viewMinSlider, 'Max', handles.dataMax);
    set (handles.viewMinSlider, 'Value', handles.dataMin);
    set (handles.viewMaxSlider, 'Max', handles.dataMax);
    set (handles.viewMaxSlider, 'Min', handles.dataMin);
    set (handles.viewMaxSlider, 'Value', handles.dataMax);  
    handles.x = round (size (handles.data, 3) / 2);
    handles.y = round (size (handles.data, 2) / 2);
    handles.z = round (size (handles.data, 1) / 2);
  updateAxes (handles);
  guidata(hObject, handles);


% --- Executes on button press in tmoviePushbutton.
function tmoviePushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to tmoviePushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  oldt = handles.t;
  for t=1:size(handles.data,4)
    handles.t = t;
    updateAxes (handles);
    uiwait (gcf, 0.5);
  end
  handles.t = oldt;
  updateAxes (handles);
  guidata(hObject, handles);


% --- Executes on button press in voiRadiobutton.
function voiRadiobutton_Callback(hObject, eventdata, handles)
% hObject    handle to voiRadiobutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of voiRadiobutton
  updateAxes (handles);
  guidata(hObject, handles);
  


% --- Executes on button press in fitRadiobutton.
function fitRadiobutton_Callback(hObject, eventdata, handles)
% hObject    handle to fitRadiobutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of fitRadiobutton
  updateAxes (handles);
  guidata(hObject, handles);


% --- Executes on selection change in voxDataPopupmenu.
function voxDataPopupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to voxDataPopupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns voxDataPopupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from voxDataPopupmenu



% --- Executes on selection change in voxErrorPopupmenu.
function voxErrorPopupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to voxErrorPopupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns voxErrorPopupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from voxErrorPopupmenu



function filterEdit_Callback(hObject, eventdata, handles)
% hObject    handle to filterEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of filterEdit as text
%        str2double(get(hObject,'String')) returns contents of filterEdit as a double


% --- Executes on selection change in filterTypePopupmenu.
function filterTypePopupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to filterTypePopupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns filterTypePopupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from filterTypePopupmenu


% --- Executes on button press in filterPushbutton.
function filterPushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to filterPushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  t = get(handles.filterTypePopupmenu,'String');
  filterType = t{get(handles.filterTypePopupmenu,'Value')}
  n = fix (str2double(get(handles.filterEdit,'String')))
  if n > 0
    save show4d_data -struct handles data;
    if strcmp (filterType, 'box') || strcmp (filterType, 'gaussian')
      for t=1:size(handles.data,4)
        handles.data(:,:,:,t) = smooth3 (squeeze (handles.data(:,:,:,t)), filterType, n);
      end
    else
      errordlg ('unknown filter', 'show4d ERROR', 'modal');
    end
  else
    errordlg ('illegal filter size', 'show4d ERROR', 'modal');
  end
  updateAxes (handles);
  guidata(hObject, handles);



% --- Executes on button press in voiLoadPushbutton.
function voiLoadPushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to voiLoadPushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  [ filename, pathname ] = uigetfile ('*.hdr');
  if ~isequal(filename, 0)
    set (handles.voiRadiobutton, 'Value', 1);
    fullname = [ pathname, filename ];
    handles.voi = voxulus ('read', 'voi', fullname);
    set (handles.msgText, 'String', sprintf ('show4d: VOI loaded from ''%s'' (%d)', filename, nnz (handles.voi)));
    updateAxes (handles);
    guidata(hObject, handles);
  end
