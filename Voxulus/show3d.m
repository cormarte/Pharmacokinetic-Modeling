function varargout = show3d(varargin)
% SHOW3D M-file for show3d.fig
%      SHOW3D, by itself, creates a new SHOW3D or raises the existing
%      singleton*.
%
%      H = SHOW3D returns the handle to a new SHOW3D or the handle to
%      the existing singleton*.
%
%      SHOW3D('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SHOW3D.M with the given input arguments.
%
%      SHOW3D('Property','Value',...) creates a new SHOW3D or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before show3d_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to show3d_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help show3d

% Last Modified by GUIDE v2.5 24-Feb-2005 09:28:11

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @show3d_OpeningFcn, ...
                   'gui_OutputFcn',  @show3d_OutputFcn, ...
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
  set (gcf, 'Name', sprintf ('show3d: z=%.0f y=%.0f x=%.0f', handles.z, handles.y, handles.x));
  set (handles.xSlider, 'Min', 1);
  nx = size(handles.data,3);
  ny = size(handles.data,2);
  nz = size(handles.data,1);
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
  set (handles.interpolateEdit, 'String', sprintf ('%.0f', max (size (handles.data))));
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
  axes (handles.colAxes);
  img = imagesc (squeeze(handles.data(handles.z,:,:)), [ viewMin viewMax ]);
  set (gca, 'Visible', 'off');
  set (img, 'Visible', 'off');
  colorbar ('WestOutside');
  axes (handles.xyAxes);
  d = squeeze (handles.data(handles.z,:,:));
  v = logical (zeros (size (d)));
  if get (handles.voiRadiobutton, 'Value')
    v = squeeze (handles.voi(handles.z,:,:));
  else
    v(handles.y,:) = 1;
    v(:,handles.x) = 1;
  end
  show_image (d, [ viewMin viewMax ], v);
  axis xy;
  set (get (handles.xyAxes, 'Children'), 'ButtonDownFcn', {@xyAxes_ButtonDownFcn,handles});
  axes (handles.yzAxes);
  d = squeeze (handles.data(:,:,handles.x));
  v = logical (zeros (size (d)));
  if get (handles.voiRadiobutton, 'Value')
    v = squeeze (handles.voi(:,:,handles.x));
  else
    v(handles.z,:) = 1;
    v(:,handles.y) = 1;
  end
  show_image (d', [ viewMin viewMax ], v');
  axis xy;
  set (get (handles.yzAxes, 'Children'), 'ButtonDownFcn', {@yzAxes_ButtonDownFcn,handles});
  axes (handles.xzAxes);
  d = squeeze (handles.data(:,handles.y,:));
  v = logical (zeros (size (d)));
  if get (handles.voiRadiobutton, 'Value')
    v = squeeze (handles.voi(:,handles.y,:));
  else
    v(handles.z,:) = 1;
    v(:,handles.x) = 1;
  end
  show_image (d, [ viewMin viewMax ], v);
  set (get (handles.xzAxes, 'Children'), 'ButtonDownFcn', {@xzAxes_ButtonDownFcn,handles});
  updateVolAxes (handles);

function updateVolAxes (handles)
  axes (handles.volAxes);
  nx = size(handles.data,3);
  ny = size(handles.data,2);
  nz = size(handles.data,1);
  viewMin = get (handles.viewMinSlider, 'Value');
  viewMax = get (handles.viewMaxSlider, 'Value');
  data = permute (handles.data, [ 2 3 1 ]);
  if strcmp (handles.volMode, 'slice')
    set (handles.volAxes, 'Visible', 'on');
    d = ones (3, 1);
    if max (size (data)) > 40
    %
    % do some smoothing and subsampling for large data sets
    %
    d = ceil (size (data) / 40);
    data = smooth3 (data, 'box', 2*fix(d/2)+1);
    data = data(1:d(1):ny,1:d(2):nx,1:d(3):nz);
    end
    slice (double (data), handles.x / d(1), handles.y / d(2), handles.z / d(3));
    axis ([ 1 size(data,1) 1 size(data,2) 1 size(data,3) viewMin viewMax ]);
    xlabel ('x');
    ylabel ('y');
    zlabel ('z');
  elseif strcmp (handles.volMode, 'rayce')
    set (handles.volAxes, 'Visible', 'on');
    set (handles.volrxSlider, 'Visible', 'on');
    set (handles.volrySlider, 'Visible', 'on');
    set (handles.volscaleSlider, 'Visible', 'on');
    rayce ('set', 'rot', [ handles.rx handles.ry handles.rz ]);
    rayce ('set', 'scale', handles.volScale);
    data = int16 (32767 * data);
    imagesc (rayce (data));
  elseif strcmp (handles.volMode, 'none')
    cla;
    set (handles.volAxes, 'Visible', 'off');
    set (handles.volrxSlider, 'Visible', 'off');
    set (handles.volrySlider, 'Visible', 'off');
    set (handles.volscaleSlider, 'Visible', 'off');
  end
  

function [ x1, x2, y1, y2 ] = getRect (ax, xmax, ymax)
  r = getrect (ax);
  x1 = round (r(1));
  x2 = round (r(1) + r(3));
  x1 = max ([ x1 1 ]);
  x2 = min ([ x2 xmax ]);
  y1 = round (r(2));
  y2 = round (r(2) + r(4));
  y1 = max ([ y1 1 ]);
  y2 = min ([ y2 ymax ]);
  

function [ x, y ] = getLine (ax, xmax, ymax)
  [ x, y ] = getline (ax);

function y = setLine (s, lx, ly, v)
  y = s;
  for i=2:size(lx,1)
    x1 = lx(i-1);
    x2 = lx(i);
    y1 = ly(i-1);
    y2 = ly(i);
    d = round (1.5 * max ([ abs(x2-x1) abs(y2-y1) ]));
    d = max ([ d 1 ]);
    dx = x2 - x1;
    dy = y2 - y1;
    for l=0:d
      a = round (x1 + l * dx / d);
      b = round (y1 + l * dy / d);
      y(b,a) = v;
    end
  end


function y = setEllipse (s, x1, x2, y1, y2, v)
  y = s;
  cx = (x1 + x2) / 2;
  cy = (y1 + y2) / 2;
  rx = (x2 - x1) / 2;
  ry = (y2 - y1) / 2;
  for w=-90:1:90
    angle = w * pi / 180;
    i2 = round (cx + rx * cos (angle));
    i1 = round (cx - rx * cos (angle));
    j = round (cy + ry * sin (angle));
    y(j,i1:i2) = v;
  end

function v = growVOI (s)
  v = logical (zeros (size (s)));
  [ y, x ] = find (s > 0);
  for i=1:size(x,1)
    v(y(i),x(i)) = 1;
    v(y(i)-1,x(i)) = 1;
    v(y(i),x(i)-1) = 1;
    v(y(i)+1,x(i)) = 1;
    v(y(i),x(i)+1) = 1;
  end

function v = rotateVOI (s)
  cx = 1 + size(s,2) / 2;
  cy = 1 + size(s,1) / 2;
  v = logical (zeros (size (s)));
  [ y, x ] = find (s == 1);
  for i=1:length(x)
    r = sqrt((x(i) - cx)^2 + (y(i) - cy)^2);
    for w=0:1:360
      angle = w * pi / 180;
      px = round (cx + r * cos (angle));
      py = round (cy + r * sin (angle));
      v(py,px) = 1;
    end
  end

% --- Executes just before show3d is made visible.
function show3d_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to show3d (see VARARGIN)

% Choose default command line output for show3d
handles.output = hObject;
handles.data = varargin{1};
if length (varargin) > 1
  handles.voi = varargin{2};
else
  handles.voi = logical (zeros (size (handles.data)));;
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
handles.x = round (nx / 2);
handles.y = round (ny / 2);
handles.z = round (nz / 2);
handles.volMode = 'none';
handles.voiMode = '+point';
handles.rx = 0;
handles.ry = 0;
handles.rz = 0;
handles.volScale = 1.0;
updateAxes (handles);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes show3d wait for user response (see UIRESUME)
%uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = show3d_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = hObject;


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


function xEditText_Callback(hObject, eventdata, handles)
% hObject    handle to xEditText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xEditText as text
%        str2double(get(hObject,'String')) returns contents of xEditText as a double


function yEditText_Callback(hObject, eventdata, handles)
% hObject    handle to yEditText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of yEditText as text
%        str2double(get(hObject,'String')) returns contents of yEditText as a double


function zEditText_Callback(hObject, eventdata, handles)
% hObject    handle to zEditText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of zEditText as text
%        str2double(get(hObject,'String')) returns contents of zEditText as a double


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


% --- Executes on button press in centerPushbutton.
function centerPushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to centerPushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  n = size (handles.data);
  [ zmax, handles.z ] = max (sum (reshape (handles.data, n(1), n(2)*n(3))'));
  [ ymax, handles.y ] = max (sum (reshape (permute (handles.data, [ 2 1 3 ]), n(2), n(1)*n(3))'));
  [ xmax, handles.x ] = max (sum (reshape (permute (handles.data, [ 3 1 2 ]), n(3), n(1)*n(2))'));
  set (handles.xSlider, 'Value', handles.x);
  set (handles.ySlider, 'Value', handles.y);
  set (handles.zSlider, 'Value', handles.z);
  updateAxes (handles);
  guidata(hObject, handles);


% --- Executes on mouse press over axes background.
function xyAxes_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to xyAxes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  p = get (gca, 'CurrentPoint');
  x = round (p(1,1));
  y = round (p(1,2));
  if get (handles.voiRadiobutton, 'Value')
    save show3d_voi -struct handles voi;
    t = get (handles.voimodePopupmenu,'String');
    handles.voiMode = t{get(handles.voimodePopupmenu,'Value')};
    if strcmp (handles.voiMode, '+point')
      handles.voi(handles.z,y,x) = 1;
    elseif strcmp (handles.voiMode, '-point')
      handles.voi(handles.z,y,x) = 0;
    elseif strcmp (handles.voiMode, '+rect')
      [ x1, x2, y1, y2 ] = getRect (gca, size(handles.data,3), size(handles.data,2));
      handles.voi(handles.z,y1:y2,x1:x2) = 1;
    elseif strcmp (handles.voiMode, '-rect')
      [ x1, x2, y1, y2 ] = getRect (gca, size(handles.data,3), size(handles.data,2));
      handles.voi(handles.z,y1:y2,x1:x2) = 0;
    elseif strcmp (handles.voiMode, '+ellipse')
      [ x1, x2, y1, y2 ] = getRect (gca, size(handles.data,3), size(handles.data,2));
      handles.voi(handles.z,:,:) = setEllipse (squeeze(handles.voi(handles.z,:,:)), x1, x2, y1, y2, 1);
    elseif strcmp (handles.voiMode, '-ellipse')
      [ x1, x2, y1, y2 ] = getRect (gca, size(handles.data,3), size(handles.data,2));
      handles.voi(handles.z,:,:) = setEllipse (squeeze(handles.voi(handles.z,:,:)), x1, x2, y1, y2, 0);
    elseif strcmp (handles.voiMode, '+line')
      [ x, y ] = getLine (gca, size(handles.data,3), size(handles.data,2));
      handles.voi(handles.z,:,:) = setLine (squeeze(handles.voi(handles.z,:,:)), x, y, 1);
    elseif strcmp (handles.voiMode, '-line')
      [ x, y ] = getLine (gca, size(handles.data,3), size(handles.data,2));
      handles.voi(handles.z,:,:) = setLine (squeeze(handles.voi(handles.z,:,:)), x, y, 0);
    elseif strcmp (handles.voiMode, '+thresh')
      t = handles.data(handles.z,y,x);
      set (handles.voithreshEdit, 'String', sprintf ('%g', t));
      k = handles.data >= t;
      handles.voi(k) = 1;
    elseif strcmp (handles.voiMode, '-thresh')
      t = handles.data(handles.z,y,x);
      set (handles.voithreshEdit, 'String', sprintf ('%g', t));
      k = handles.data <= t;
      handles.voi = 1;
    elseif strcmp (handles.voiMode, 'grow')
      handles.voi(handles.z,:,:) = growVOI (squeeze(handles.voi(handles.z,:,:)));
    elseif strcmp (handles.voiMode, 'rotate')
      handles.voi(handles.z,:,:) = rotateVOI (squeeze(handles.voi(handles.z,:,:)));
    elseif strcmp (handles.voiMode, 'toggle')
      handles.voi(handles.z,y,x) = 1 - handles.voi(handles.z,y,x);
    elseif strcmp (handles.voiMode, 'clear')
      handles.voi(handles.z,:,:) = 0;
    end
  else
    handles.x = x;
    handles.y = y;
    set (handles.xSlider, 'Value', handles.x);
    set (handles.ySlider, 'Value', handles.y);
  end
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
  if get (handles.voiRadiobutton, 'Value')
    save show3d_voi -struct handles voi;
    t = get (handles.voimodePopupmenu,'String');
    handles.voiMode = t{get(handles.voimodePopupmenu,'Value')};
    if strcmp (handles.voiMode, '+point')
      handles.voi(z,y,handles.x) = 1;
    elseif strcmp (handles.voiMode, '-point')
      handles.voi(z,y,handles.x) = 0;
    elseif strcmp (handles.voiMode, '+rect')
      [ z1, z2, y1, y2 ] = getRect (gca, size(handles.data,1), size(handles.data,2));
      handles.voi(z1:z2,y1:y2,handles.x) = 1;
    elseif strcmp (handles.voiMode, '-rect')
      [ z1, z2, y1, y2 ] = getRect (gca, size(handles.data,1), size(handles.data,2));
      handles.voi(z1:z2,y1:y2,handles.x) = 0;
    elseif strcmp (handles.voiMode, '+ellipse')
      [ z1, z2, y1, y2 ] = getRect (gca, size(handles.data,1), size(handles.data,2));
      handles.voi(:,:,handles.x) = setEllipse (squeeze(handles.voi(:,:,handles.x)), y1, y2, z1, z2, 1);
    elseif strcmp (handles.voiMode, '-ellipse')
      [ z1, z2, y1, y2 ] = getRect (gca, size(handles.data,1), size(handles.data,2));
      handles.voi(:,:,handles.x) = setEllipse (squeeze(handles.voi(:,:,handles.x)), y1, y2, z1, z2, 0);
    elseif strcmp (handles.voiMode, '+line')
      [ z, y ] = getLine (gca, size(handles.data,1), size(handles.data,2));
      handles.voi(:,:,handles.x) = setLine (squeeze(handles.voi(:,:,handles.x)), y, z, 1);
    elseif strcmp (handles.voiMode, '-line')
      [ z, y ] = getLine (gca, size(handles.data,1), size(handles.data,2));
      handles.voi(:,:,handles.x) = setLine (squeeze(handles.voi(:,:,handles.x)), y, z, 0);
    elseif strcmp (handles.voiMode, '+thresh')
      t = handles.data(z,y,handles.x);
      set (handles.voithreshEdit, 'String', sprintf ('%g', t));
      k = handles.data >= t;
      handles.voi(k) = 1;
    elseif strcmp (handles.voiMode, '-thresh')
      t = handles.data(z,y,handles.x);
       set (handles.voithreshEdit, 'String', sprintf ('%g', t));
      k = handles.data <= t;
      handles.voi(k) = 1;
    elseif strcmp (handles.voiMode, 'grow')
      handles.voi(:,:,handles.x) = growVOI (squeeze(handles.voi(:,:,handles.x)));
    elseif strcmp (handles.voiMode, 'rotate')
      handles.voi(:,:,handles.x) = rotateVOI (squeeze(handles.voi(:,:,handles.x)));
    elseif strcmp (handles.voiMode, 'toggle')
      handles.voi(z,y,handles.x) = 1 - handles.voi(z,y,handles.x);
    elseif strcmp (handles.voiMode, 'clear')
      handles.voi(:,:,handles.x) = 0;
    end
  else
    handles.z = z;
    handles.y = y;
    set (handles.zSlider, 'Value', handles.z);
    set (handles.ySlider, 'Value', handles.y);
  end
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
  if get (handles.voiRadiobutton, 'Value')
    save show3d_voi -struct handles voi;
    t = get (handles.voimodePopupmenu,'String');
    handles.voiMode = t{get(handles.voimodePopupmenu,'Value')};
    if strcmp (handles.voiMode, '+point')
      handles.voi(z,handles.y,x) = 1;
    elseif strcmp (handles.voiMode, '-point')
      handles.voi(z,handles.y,x) = 0;
    elseif strcmp (handles.voiMode, '+rect')
      [ x1, x2, z1, z2 ] = getRect (gca, size(handles.data,3), size(handles.data,1));
      handles.voi(z1:z2,handles.y,x1:x2) = 1;
    elseif strcmp (handles.voiMode, '-rect')
      [ x1, x2, z1, z2 ] = getRect (gca, size(handles.data,3), size(handles.data,1));
      handles.voi(z1:z2,handles.y,x1:x2) = 0;
    elseif strcmp (handles.voiMode, '+ellipse')
      [ x1, x2, y1, y2 ] = getRect (gca, size(handles.data,3), size(handles.data,1));
      handles.voi(:,handles.y,:) = setEllipse (squeeze(handles.voi(:,handles.x,:)), x1, x2, y1, y2, 1);
    elseif strcmp (handles.voiMode, '-ellipse')
      [ x1, x2, y1, y2 ] = getRect (gca, size(handles.data,3), size(handles.data,1));
      handles.voi(:,handles.y,:) = setEllipse (squeeze(handles.voi(:,handles.y,:)), x1, x2, y1, y2, 0);
    elseif strcmp (handles.voiMode, '+line')
      [ x, z ] = getLine (gca, size(handles.data,3), size(handles.data,1));
      handles.voi(:,handles.y,:) = setLine (squeeze(handles.voi(:,handles.y,:)), x, z, 1);
    elseif strcmp (handles.voiMode, '-line')
      [ x, z ] = getLine (gca, size(handles.data,3), size(handles.data,1));
      handles.voi(:,handles.y,:) = setLine (squeeze(handles.voi(:,handles.y,:)), x, z, 0);
    elseif strcmp (handles.voiMode, '+thresh')
      t = handles.data(z,handles.y,x);
      set (handles.voithreshEdit, 'String', sprintf ('%g', t));
      k = handles.data >= t;
      handles.voi(k) = 1;
    elseif strcmp (handles.voiMode, '-thresh')
      t = handles.data(z,handles.y,x);
       set (handles.voithreshEdit, 'String', sprintf ('%g', t));
      k = handles.data <= t;
      handles.voi(k) = 1;
    elseif strcmp (handles.voiMode, 'grow')
      handles.voi(:,handles.y,:) = growVOI (squeeze(handles.voi(:,handles.y,:)));
    elseif strcmp (handles.voiMode, 'rotate')
      handles.voi(:,handles.y,:) = rotateVOI (squeeze(handles.voi(:,handles.y,:)));
    elseif strcmp (handles.voiMode, 'toggle')
      handles.voi(z,handles.y,x) = 1 - handles.voi(z,handles.y,x);
    elseif strcmp (handles.voiMode, 'clear')
      handles.voi(:,handles.y,:) = 0;
    end
  else
    handles.x = x;
    handles.z = z;
    set (handles.xSlider, 'Value', handles.x);
    set (handles.zSlider, 'Value', handles.z);
  end
  updateAxes (handles);
  guidata (gcf, handles);


% --- Executes on button press in rotateZ_Pushbutton.
function rotateZ_Pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to rotateZ_Pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  rz = str2double (get (handles.rzEdit, 'String'));
  if rz == 0
    [ x, y ] = getline (handles.xyAxes);
    if length (x) == 2 && length (y) == 2
      dx = x(2) - x(1);
      dy = y(2) - y(1);
      rz = 90.0 - atan2 (dy, dx) * 180.0 / pi;
    end
  end
  if rz ~= 0
    save show3d_data -struct handles data;
    save show3d_voi -struct handles voi;
    % rotation algorithm requires cubic volume: perform zero padding if
    % necessary
    handles = zeroPadding ( handles );
    c = [ handles.x, handles.y, handles.z ];
    fprintf ('''rotate'', [ 0.0 0.0 %.1f ], [ %.1f %.1f %.1f ]\n', rz, c);
    handles.data = permute (rotate3 (permute (handles.data, [ 3 2 1 ]), 0.0, 0.0, rz, c), [ 3 2 1]);
    set (handles.rzEdit, 'String', '0');
  end
  updateAxes (handles);
  guidata(hObject, handles);


% --- Executes on button press in rotateX_Pushbutton.
function rotateX_Pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to rotateX_Pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  rx = str2double (get (handles.rxEdit, 'String'));
  if rx == 0
    [ x, y ] = getline (handles.yzAxes);
    if length (x) == 2 && length (y) == 2
      dx = x(2) - x(1);
      dy = y(2) - y(1);
      rx = 90.0 - atan2 (dy, dx) * 180.0 / pi;
    end
  end
  if rx ~= 0
    save show3d_data -struct handles data;
    save show3d_voi -struct handles voi;
    % rotation algorithm requires cubic volume: perform zero padding if
    % necessary
    handles = zeroPadding ( handles );
    c = [ handles.x, handles.y, handles.z ];
    fprintf ('''rotate'', [ 0.0 %.1f 0.0 ], [ %.1f %.1f %.1f ]\n', rx, c);
    handles.data = permute (rotate3 (permute (handles.data, [ 3 2 1 ]), 0.0, rx, 0.0, c), [ 3 2 1 ]);
    set (handles.rxEdit, 'String', '0');
  end
  updateAxes (handles);
  guidata(hObject, handles);


% --- Executes on button press in rotateY_Pushbutton.
function rotateY_Pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to rotateY_Pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  ry = str2double (get (handles.ryEdit, 'String'));
  if ry == 0
    [ x, y ] = getline (handles.xzAxes);
    if length (x) == 2 && length (y) == 2
      dx = x(2) - x(1);
      dy = y(2) - y(1);
      ry = atan2 (dy, dx) * 180.0 / pi - 90.0;
    end
  end
  if ry ~= 0
    save show3d_data -struct handles data;
    save show3d_voi -struct handles voi;
    % rotation algorithm requires cubic volume: perform zero padding if
    % necessary
    handles = zeroPadding ( handles );
    c = [ handles.x, handles.y, handles.z ];
    fprintf ('''rotate'', [ %.1f 0.0 0.0 ], [ %.1f %.1f %.1f ]\n', ry, c);
    handles.data = permute (rotate3 (permute (handles.data, [ 3 2 1 ]), ry, 0.0, 0.0, c), [ 3 2 1 ]);
    set (handles.ryEdit, 'String', '0');
  end
  updateAxes (handles);
  guidata(hObject, handles);


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


% --- Executes during object creation, after setting all properties.
function viewMinSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to viewMinSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background, change
%       'usewhitebg' to 0 to use default.  See ISPC and COMPUTER.
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes during object creation, after setting all properties.
function viewMaxSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to viewMaxSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background, change
%       'usewhitebg' to 0 to use default.  See ISPC and COMPUTER.
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on button press in cropxyPushbutton.
function cropxyPushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to cropxyPushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  rxy = getrect (handles.xyAxes);
  x1 = round (rxy(1));
  x2 = round (rxy(1) + rxy(3));
  x1 = max ([ x1 1 ]);
  x2 = min ([ x2 size(handles.data,3) ]);
  y1 = round (rxy(2));
  y2 = round (rxy(2) + rxy(4));
  y1 = max ([ y1 1 ]);
  y2 = min ([ y2 size(handles.data,2) ]);
  if (rxy(3) > 1) && (rxy(4) > 1)
    save show3d_data -struct handles data
    save show3d_voi -struct handles voi;
    set (handles.msgText, 'String', sprintf ('show3d: cropping from %dx%dx%d to %dx%dx%d', ...
      size (handles.data, 3), size (handles.data, 2), size (handles.data, 1), ...
      x2 - x1 + 1, y2 - y1 + 1, size (handles.data, 1)));
    z1 = 1;
    z2 = size (handles.data, 1);
    fprintf ('''crop'', %d:%d, %d:%d, %d:%d\n', z1, z2, y1, y2, x1, x2);
    handles.data = handles.data(:,y1:y2,x1:x2);
    handles.voi = handles.voi(:,y1:y2,x1:x2);
    handles.dataMin = min (handles.data(:));
    handles.dataMax = max (handles.data(:));
    set (handles.viewMinSlider, 'Min', handles.dataMin);
    set (handles.viewMinSlider, 'Max', handles.dataMax);
    set (handles.viewMinSlider, 'Value', handles.dataMin);
    set (handles.viewMaxSlider, 'Max', handles.dataMax);
    set (handles.viewMaxSlider, 'Min', handles.dataMin);
    set (handles.viewMaxSlider, 'Value', handles.dataMax);
    handles.x = handles.x - x1 + 1;
    handles.y = handles.y - y1 + 1;
    handles.z = handles.z - z1 + 1;
    if ( ( handles.x < 1 ) || ( handles.x > x2 - x1 + 1 ) )
        handles.x = round ( ( x2 - x1 + 1 ) / 2 );
    end 
    if ( ( handles.y < 1 ) || ( handles.y > y2 - y1 + 1 ) )
        handles.y = round ( ( y2 - y1 + 1 ) / 2 );
    end 
    if ( ( handles.z < 1 ) || ( handles.z > z2 - z1 + 1 ) )
        handles.z = round ( ( z2 - z1 + 1 ) / 2 );
    end 
    updateAxes (handles);
    guidata(hObject, handles);
  end


% --- Executes on button press in rotate3D_Radiobutton.
function rotate3D_Radiobutton_Callback(hObject, eventdata, handles)
% hObject    handle to rotate3D_Radiobutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rotate3D_Radiobutton
  axes (handles.volAxes);
  if get(hObject, 'Value')
    rotate3d on
  else
    rotate3d off
  end


function interpolateEdit_Callback(hObject, eventdata, handles)
% hObject    handle to interpolateEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of interpolateEdit as text
%        str2double(get(hObject,'String')) returns contents of interpolateEdit as a double


% --- Executes on button press in interpolatePushbutton.
function interpolatePushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to interpolatePushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  save show3d_data -struct handles data
  n = str2double (get (handles.interpolateEdit, 'String'));
  nx = size(handles.data,3);
  ny = size(handles.data,2);
  nz = size(handles.data,1);
  [ x, y, z ] = ndgrid (1:nx, 1:ny, 1:nz);
  [ xi, yi, zi ] = meshgrid (linspace (1, nx, n), linspace (1, ny, n), linspace (1, nz, n));
  fprintf ('''interp'', [ %d %d %d ]\n', n, n, n);
  tmp = interp3 (x, y, z, permute (handles.data, [ 3 2 1 ]), xi, yi, zi, 'cubic', min (handles.data(:)));
  handles.data = permute (tmp, [ 3 2 1 ]);
  handles.voi = logical (zeros (size (handles.data)));;
  handles.x = round (n / 2);
  handles.y = round (n / 2);
  handles.z = round (n / 2);
  updateAxes (handles);
  guidata (hObject, handles);


% --- Executes on button press in cropyzPushbutton.
function cropyzPushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to cropyzPushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  r = getrect (handles.yzAxes);
  z1 = round (r(1));
  z2 = round (r(1) + r(3));
  z1 = max ([ z1 1 ]);
  z2 = min ([ z2 size(handles.data,1) ]);
  y1 = round (r(2));
  y2 = round (r(2) + r(4));
  y1 = max ([ y1 1 ]);
  y2 = min ([ y2 size(handles.data,2) ]);
  if (r(3) > 1) && (r(4) > 1)
    save show3d_data -struct handles data;
    save show3d_voi -struct handles voi;
    set (handles.msgText, 'String', sprintf ('show3d: cropping from %dx%dx%d to %dx%dx%d', ...
      size (handles.data, 3), size (handles.data, 2), size (handles.data, 1), ...
      size (handles.data, 3), y2 - y1 + 1, z2 - z1 + 1));
    x1 = 1;
    x2 = size (handles.data, 3);
    fprintf ('''crop'', %d:%d, %d:%d, %d:%d\n', z1, z2, y1, y2, x1, x2);
    handles.data = handles.data(z1:z2,y1:y2,:);
    handles.voi = handles.voi(z1:z2,y1:y2,:);
    handles.dataMin = min (handles.data(:));
    handles.dataMax = max (handles.data(:));
    set (handles.viewMinSlider, 'Min', handles.dataMin);
    set (handles.viewMinSlider, 'Max', handles.dataMax);
    set (handles.viewMinSlider, 'Value', handles.dataMin);
    set (handles.viewMaxSlider, 'Max', handles.dataMax);
    set (handles.viewMaxSlider, 'Min', handles.dataMin);
    set (handles.viewMaxSlider, 'Value', handles.dataMax);
    handles.x = handles.x - x1 + 1;
    handles.y = handles.y - y1 + 1;
    handles.z = handles.z - z1 + 1;
    if ( ( handles.x < 1 ) || ( handles.x > x2 - x1 + 1 ) )
        handles.x = round ( ( x2 - x1 + 1 ) / 2 );
    end 
    if ( ( handles.y < 1 ) || ( handles.y > y2 - y1 + 1 ) )
        handles.y = round ( ( y2 - y1 + 1 ) / 2 );
    end 
    if ( ( handles.z < 1 ) || ( handles.z > z2 - z1 + 1 ) )
        handles.z = round ( ( z2 - z1 + 1 ) / 2 );
    end 
    updateAxes (handles);
    guidata(hObject, handles);
  end

% --- Executes on button press in cropxzPushbutton.
function cropxzPushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to cropxzPushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  r = getrect (handles.xzAxes);
  x1 = round (r(1));
  x2 = round (r(1) + r(3));
  x1 = max ([ x1 1 ]);
  x2 = min ([ x2 size(handles.data,3) ]);
  z1 = round (r(2));
  z2 = round (r(2) + r(4));
  z1 = max ([ z1 1 ]);
  z2 = min ([ z2 size(handles.data,1) ]);
  if (r(3) > 1) && (r(4) > 1)
    save show3d_data -struct handles data
    save show3d_voi -struct handles voi;
    set (handles.msgText, 'String', sprintf ('show3d: cropping from %dx%dx%d to %dx%dx%d', ...
      size (handles.data, 3), size (handles.data, 2), size (handles.data, 1), ...
      x2 - x1 + 1, size (handles.data, 2), z2 - z1 + 1));
    y1 = 1;
    y2 = size (handles.data, 2);
    fprintf ('''crop'', %d:%d, %d:%d, %d:%d\n', z1, z2, y1, y2, x1, x2);
    handles.data = handles.data(z1:z2,:,x1:x2);
    handles.voi = handles.voi(z1:z2,:,x1:x2);
    handles.dataMin = min (handles.data(:));
    handles.dataMax = max (handles.data(:));
    set (handles.viewMinSlider, 'Min', handles.dataMin);
    set (handles.viewMinSlider, 'Max', handles.dataMax);
    set (handles.viewMinSlider, 'Value', handles.dataMin);
    set (handles.viewMaxSlider, 'Max', handles.dataMax);
    set (handles.viewMaxSlider, 'Min', handles.dataMin);
    set (handles.viewMaxSlider, 'Value', handles.dataMax);  
    handles.x = handles.x - x1 + 1;
    handles.y = handles.y - y1 + 1;
    handles.z = handles.z - z1 + 1;
    if ( ( handles.x < 1 ) || ( handles.x > x2 - x1 + 1 ) )
        handles.x = round ( ( x2 - x1 + 1 ) / 2 );
    end 
    if ( ( handles.y < 1 ) || ( handles.y > y2 - y1 + 1 ) )
        handles.y = round ( ( y2 - y1 + 1 ) / 2 );
    end 
    if ( ( handles.z < 1 ) || ( handles.z > z2 - z1 + 1 ) )
        handles.z = round ( ( z2 - z1 + 1 ) / 2 );
    end 
    updateAxes (handles);
    guidata(hObject, handles);
  end

% --- Executes on button press in flipxyPushbutton.
function flipxyPushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to flipxyPushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  save show3d_data -struct handles data;
  save show3d_voi -struct handles voi;
  handles.data = permute (handles.data, [ 1 3 2 ]);
  handles.voi = permute (handles.voi, [ 1 3 2 ]);
  fprintf ('''flipxy''\n');
  handles.x = round (size (handles.data, 3) / 2);
  handles.y = round (size (handles.data, 2) / 2);
  handles.z = round (size (handles.data, 1) / 2);
  updateAxes (handles);
  guidata(hObject, handles);

% --- Executes on button press in flipxzPushbutton.
function flipxzPushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to flipxzPushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  save show3d_data -struct handles data;
  save show3d_voi -struct handles voi;
  handles.data = permute (handles.data, [ 3 2 1 ]);
  handles.voi = permute (handles.voi, [ 3 2 1 ]);
  fprintf ('''flipxz''\n');
  handles.x = round (size (handles.data, 3) / 2);
  handles.y = round (size (handles.data, 2) / 2);
  handles.z = round (size (handles.data, 1) / 2);
  updateAxes (handles);
  guidata(hObject, handles);


% --- Executes on button press in flipyzPushbutton.
function flipyzPushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to flipyzPushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  save show3d_data -struct handles data;
  save show3d_voi -struct handles voi;
  handles.data = permute (handles.data, [ 2 3 1 ]);
  handles.voi = permute (handles.voi, [ 2 3 1 ]);
  fprintf ('''flipyz''\n');
  handles.x = round (size (handles.data, 3) / 2);
  handles.y = round (size (handles.data, 2) / 2);
  handles.z = round (size (handles.data, 1) / 2);
  updateAxes (handles);
  guidata(hObject, handles);


% --- Executes on button press in revxPushbutton.
function revxPushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to revxPushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  save show3d_data -struct handles data;
  save show3d_voi -struct handles voi;
  nx = size(handles.data,3);
  for z = 1:size(handles.data,1)
    for y = 1:size(handles.data,2)
      handles.data(z,y,:) = handles.data(z,y,nx:-1:1);
      handles.voi(z,y,:) = handles.voi(z,y,nx:-1:1);
    end;
  end
  fprintf ('''revx''\n');
  updateAxes (handles);
  guidata(hObject, handles);


% --- Executes on button press in revyPushbutton.
function revyPushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to revyPushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  save show3d_data -struct handles data;
  save show3d_voi -struct handles voi;
  ny = size(handles.data,2);
  for z = 1:size(handles.data,1)
    for x = 1:size(handles.data,3)
      handles.data(z,:,x) = handles.data(z,ny:-1:1,x);
      handles.voi(z,:,x) = handles.voi(z,ny:-1:1,x);
    end;
  end
  fprintf ('''revy''\n');
  updateAxes (handles);
  guidata(hObject, handles);


% --- Executes on button press in revzPushbutton.
function revzPushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to revzPushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  save show3d_data -struct handles data;
  save show3d_voi -struct handles voi;
  nz = size(handles.data,1);
  for y = 1:size(handles.data,2)
    for x = 1:size(handles.data,3)
      handles.data(:,y,x) = handles.data(nz:-1:1,y,x);
      handles.voi(:,y,x) = handles.voi(nz:-1:1,y,x);
    end;
  end
  fprintf ('''revz''\n');
  updateAxes (handles);
  guidata(hObject, handles);


% --- Executes on selection change in volmodePopupmenu.
function volmodePopupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to volmodePopupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns volmodePopupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from volmodePopupmenu
  t = get(hObject,'String');
  handles.volMode = t{get(hObject,'Value')};
  updateVolAxes (handles);
  guidata(hObject, handles);

% --- Executes on slider movement.
function volrxSlider_Callback(hObject, eventdata, handles)
% hObject    handle to volrxSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
  handles.rx = get (hObject,'Value');
  updateVolAxes (handles);
  guidata(hObject, handles);


% --- Executes on button press in undoPushbutton.
function undoPushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to undoPushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  load show3d_data data;
  handles.data = data;
  load show3d_voi voi;
  handles.voi = voi;
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


% --- Executes on button press in zmoviePushButton.
function zmoviePushButton_Callback(hObject, eventdata, handles)
% hObject    handle to zmoviePushButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  oldz = handles.z;
  for z=1:size(handles.data,1)
    handles.z = z;
    updateAxes (handles);
    uiwait (gcf, 0.5);
  end
  handles.z = oldz;
  updateAxes (handles);
  guidata(hObject, handles);


% --- Executes on button press in xmoviePushbutton.
function xmoviePushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to xmoviePushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  oldx = handles.x;
  for x=1:size(handles.data,3)
    handles.x = x;
    updateAxes (handles);
    uiwait (gcf, 0.5);
  end
  handles.x = oldx;
  updateAxes (handles);
  guidata(hObject, handles);


% --- Executes on button press in ymoviePushbutton.
function ymoviePushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to ymoviePushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  oldy = handles.y;
  for y=1:size(handles.data,2)
    handles.y = y;
    updateAxes (handles);
    uiwait (gcf, 0.5);
  end
  handles.y = oldy;
  updateAxes (handles);
  guidata(hObject, handles);



function rxEdit_Callback(hObject, eventdata, handles)
% hObject    handle to rxEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rxEdit as text
%        str2double(get(hObject,'String')) returns contents of rxEdit as a double


function rzEdit_Callback(hObject, eventdata, handles)
% hObject    handle to rzEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rzEdit as text
%        str2double(get(hObject,'String')) returns contents of rzEdit as a double



function ryEdit_Callback(hObject, eventdata, handles)
% hObject    handle to ryEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ryEdit as text
%        str2double(get(hObject,'String')) returns contents of ryEdit as a double


% --- Executes on button press in voiRadiobutton.
function voiRadiobutton_Callback(hObject, eventdata, handles)
% hObject    handle to voiRadiobutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of voiRadiobutton
  if get (hObject,'Value')
    set (handles.msgText, 'String', 'show3d: VOI mode on');
  else
    set (handles.msgText, 'String', 'show3d: VOI mode off');
  end
  updateAxes (handles);
  guidata(hObject, handles);


% --- Executes on button press in voisavePushbutton.
function voisavePushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to voisavePushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  [ filename, pathname ] = uiputfile ('*.hdr');
  if ~isequal(filename, 0)
    voxulus ('write', 'voi', sprintf ('%s/%s', pathname, filename), handles.voi);
    set (handles.msgText, 'String', sprintf ('show3d: VOI saved to ''%s''', filename));
  end


% --- Executes on button press in voiloadPushbutton.
function voiloadPushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to voiloadPushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  [ filename, pathname ] = uigetfile ('*.hdr');
  if ~isequal(filename, 0)
    save show3d_voi -struct handles voi;
    set (handles.voiRadiobutton, 'Value', 1);
    fullname = [ pathname, filename ];
    handles.voi = voxulus ('read', 'voi', fullname);
    set (handles.msgText, 'String', sprintf ('show3d: VOI loaded from ''%s'' (%d)', filename, nnz (handles.voi)));
    if get (handles.voiRadiobutton, 'Value')
      updateAxes (handles);
      guidata(hObject, handles);
    end
  end


% --- Executes on button press in voixminusPushbutton.
function voixminusPushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to voixminusPushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  if get (handles.voiRadiobutton, 'Value')
    if handles.x > 1
      save show3d_voi -struct handles voi;
      handles.voi(:,:,handles.x) = handles.voi(:,:,handles.x-1); 
    end
    updateAxes (handles);
    guidata(hObject, handles);
  end


% --- Executes on button press in voixplusPushbutton.
function voixplusPushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to voixplusPushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  if get (handles.voiRadiobutton, 'Value')
    if handles.x < size(handles.data,3)
      save show3d_voi -struct handles voi;
      handles.voi(:,:,handles.x) = handles.voi(:,:,handles.x+1); 
    end
    updateAxes (handles);
    guidata(hObject, handles);
  end


% --- Executes on button press in voiyminusPushbutton.
function voiyminusPushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to voiyminusPushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  if get (handles.voiRadiobutton, 'Value')
    if handles.y > 1
      save show3d_voi -struct handles voi;
      handles.voi(:,handles.y,:) = handles.voi(:,handles.y-1,:); 
    end
    updateAxes (handles);
    guidata(hObject, handles);
  end


% --- Executes on button press in voiyplusPushbutton.
function voiyplusPushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to voiyplusPushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  if get (handles.voiRadiobutton, 'Value')
    if handles.y < size(handles.data,2)
      save show3d_voi -struct handles voi;
      handles.voi(:,handles.y,:) = handles.voi(:,handles.y+1,:); 
    end
    updateAxes (handles);
    guidata(hObject, handles);
  end


% --- Executes on button press in voizminusPushbutton.
function voizminusPushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to voizminusPushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  if get (handles.voiRadiobutton, 'Value')
    if handles.z > 1
      save show3d_voi -struct handles voi;
      handles.voi(handles.z,:,:) = handles.voi(handles.z-1,:,:); 
    end
    updateAxes (handles);
    guidata(hObject, handles);
  end


% --- Executes on button press in voizplusPushbutton.
function voizplusPushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to voizplusPushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  if get (handles.voiRadiobutton, 'Value')
    if handles.z < size(handles.data,1)
      save show3d_voi -struct handles voi;
      handles.voi(handles.z,:,:) = handles.voi(handles.z+1,:,:); 
    end
    updateAxes (handles);
    guidata(hObject, handles);
  end


% --- Executes on button press in clearvoiPushbutton.
function clearvoiPushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to clearvoiPushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  if get (handles.voiRadiobutton, 'Value')
    save show3d_voi -struct handles voi;
    handles.voi = logical (zeros (size (handles.data))); 
    updateAxes (handles);
    guidata(hObject, handles);
  end


% --- Executes on selection change in voimodePopupmenu.
function voimodePopupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to voimodePopupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns voimodePopupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from voimodePopupmenu
  t = get(hObject,'String');
  handles.voiMode = t{get(hObject,'Value')};
  set (handles.msgText, 'String', sprintf ('show3d: VOI mode ''%s''', handles.voiMode));
  guidata(hObject, handles);



function voithreshEdit_Callback(hObject, eventdata, handles)
% hObject    handle to voithreshEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of voithreshEdit as text
%        str2double(get(hObject,'String')) returns contents of voithreshEdit as a double
  save show3d_voi -struct handles voi;
  t = str2double(get(hObject,'String'));
  if t >= 0
    k = handles.data >= t;
  else
    k = handles.data <= t;
  end
  handles.voi(k) = 1;
  updateAxes (handles);
  guidata(hObject, handles);


% --- Executes on button press in voismoothPushbutton.
function voismoothPushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to voismoothPushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  save show3d_voi -struct handles voi;
  nBefore = nnz (handles.voi);
  t = str2double (get (handles.voithreshEdit, 'String'));
  if t == 0
    tmp = smooth3 (handles.voi, 'gaussian');
  else
    tmp = smooth3 (handles.voi, 'gaussian', [ 3 3 3 ], t);
  end
  handles.voi = logical (zeros (size (handles.voi)));
  handles.voi(tmp>=0.25) = 1;
  set (handles.msgText, 'String', sprintf ('show3d: VOI smoothed (%d -> %d)', nBefore, nnz (handles.voi)));
  updateAxes (handles);
  guidata(hObject, handles);


% --- Executes on button press in saveimagePushbutton.
function saveimagePushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to saveimagePushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  [ filename, pathname ] = uiputfile ('*.png');
  if ~isequal(filename, 0)
    n = round (str2double (get (handles.interpolateEdit, 'String')));
    filename = filename(1:length(filename)-4);
    axes (handles.xyAxes);
    xyImg = sprintf ('%s/%s_xy.png', pathname, filename);
    d = flipdim (get (get (gca, 'Children'), 'CData'), 1);
    imwrite (interpolate_image (d, n), xyImg);
    axes (handles.yzAxes);
    yzImg = sprintf ('%s/%s_yz.png', pathname, filename);
    d = flipdim (get (get (gca, 'Children'), 'CData'), 1);
    imwrite (interpolate_image (d, n), yzImg);
    axes (handles.xzAxes);
    xzImg = sprintf ('%s/%s_xz.png', pathname, filename);
    d = get (get (gca, 'Children'), 'CData');
    imwrite (interpolate_image (d, n), xzImg);
  end


% --- Executes on slider movement.
function volrySlider_Callback(hObject, eventdata, handles)
% hObject    handle to volrySlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
  handles.ry = get (hObject,'Value');
  updateVolAxes (handles);
  guidata(hObject, handles);


% --- Executes on slider movement.
function volscaleSlider_Callback(hObject, eventdata, handles)
% hObject    handle to volscaleSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
  handles.volScale = get (hObject,'Value');
  updateVolAxes (handles);
  guidata(hObject, handles);


% --- Executes on button press in padPushbutton.
function padPushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to padPushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  save show3d_data -struct handles data
  save show3d_voi -struct handles voi;
  % perform zero padding to achieve cubic volume
  handles = zeroPadding ( handles );
  % reflect changes on panel
  updateAxes (handles);
  guidata (hObject, handles);


function RetHandles = zeroPadding ( handles )
    [nz ny nx] = size ( handles.data );
    if ( ( nx == ny ) && ( nx == nz ) )
        % nothing to do for us
        RetHandles = handles;
        return;
    end
    % get maximum dimension
    n = max ( [nz ny nx] );
    % define starting point of original data in extended volume
    ox = fix ( 1 + (n - nx) / 2 );
    oy = fix ( 1 + (n - ny) / 2 );
    oz = fix ( 1 + (n - nz) / 2 );
    fprintf ( 1, '''pad'', %d\n', n);
    % create the new volume for data
    Cube = zeros (n, n, n);
    Cube(oz:oz+nz-1,oy:oy+ny-1,ox:ox+nx-1) = handles.data;
    handles.data = Cube;
    % create the new volume for VOI
    Cube = zeros (n, n, n);
    Cube(oz:oz+nz-1,oy:oy+ny-1,ox:ox+nx-1) = handles.voi;
    handles.voi = Cube;
    % reposition the center point according to new origin
    handles.x = handles.x + ox;
    handles.y = handles.y + oy;
    handles.z = handles.z + oz;
    % set return value before terminating!
    RetHandles = handles;

