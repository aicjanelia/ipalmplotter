function varargout = ipalmcluster(varargin)
% IPALMCLUSTER MATLAB code for ipalmcluster.fig
%      IPALMCLUSTER, by itself, creates a new IPALMCLUSTER or raises the existing
%      singleton*.
%
%      H = IPALMCLUSTER returns the handle to a new IPALMCLUSTER or the handle to
%      the existing singleton*.
%
%      IPALMCLUSTER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IPALMCLUSTER.M with the given input arguments.
%
%      IPALMCLUSTER('Property','Value',...) creates a new IPALMCLUSTER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ipalmcluster_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ipalmcluster_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ipalmcluster

% Last Modified by GUIDE v2.5 02-Jul-2019 11:29:46

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ipalmcluster_OpeningFcn, ...
                   'gui_OutputFcn',  @ipalmcluster_OutputFcn, ...
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


% --- Executes just before ipalmcluster is made visible.
function ipalmcluster_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ipalmcluster (see VARARGIN)

% Choose default command line output for ipalmcluster
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ipalmcluster wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ipalmcluster_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in loaddata.
function loaddata_Callback(hObject, eventdata, handles)
% hObject    handle to loaddata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
output = inputdlg('Pixel Size (nm)','Specify Pixel Size',1,{'133'});
pixelsize = str2double(output{1});
[data,paramlist,filename,pathname] = ipalmimport(pixelsize);
%paramlist = strsplit(paramlist,'\t');
%paramlist = paramlist(1:(end-1));
xposidx = find(contains(paramlist,'X_Position'),1);
yposidx = find(contains(paramlist,'Y_Position'),1);
zposidx = find(contains(paramlist,'Unwrapped_Z'),1);
labelidx = find(contains(paramlist,'Label_Set'),1);
clustercol = find(contains(paramlist,'Cluster_Index'),1);
xpos = data(:,xposidx);
ypos = data(:,yposidx);
zpos = data(:,zposidx);
label = data(:,labelidx);
clusteridx = data(:,clustercol);
set(handles.currfile,'String',fullfile(pathname,filename));
data = [xpos ypos zpos label clusteridx];
setappdata(handles.loaddata,'data',data);
i = 0;
if ismember(1,label)
    i = i+1;
    diststring{i} = 'Label 1 Peaks';
    clus1 = data(label == 1 & clusteridx > 0,:);
    if size(clus1,1)>0
        i = i+1;
        diststring{i} = 'Label 1 Clusters';
    end
end

if ismember(2,label)
    i = i+1;
    diststring{i} = 'Label 2 Peaks';
    clus2 = data(label == 2 & clusteridx > 0,:);
    if size(clus2,1)>0
        i = i+1;
        diststring{i} = 'Label 2 Clusters';
    end
end

if ismember(3,label)
    i = i+1;
    diststring{i} = 'Label 3 Peaks';
    clus3 = data(label == 3 & clusteridx > 0,:);
    if size(clus3,1)>0
        i = i+1;
        diststring{i} = 'Label 3 Clusters';
    end
end

handles.distfrom.String = diststring;
handles.distto.String = diststring;   
    

function currfile_Callback(hObject, eventdata, handles)
% hObject    handle to currfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of currfile as text
%        str2double(get(hObject,'String')) returns contents of currfile as a double


% --- Executes during object creation, after setting all properties.
function currfile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to currfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in plotdata.
function plotdata_Callback(hObject, eventdata, handles)
% hObject    handle to plotdata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data = getappdata(handles.loaddata,'data');
xpos = data(:,1);
ypos = data(:,2);
zpos = data(:,3);
label = data(:,4);
labels = unique(label);
colors = ['r','g','b'];
clusteridx = data(:,5);
plottype = get(handles.plottype,'Value');
axes(handles.axes1);
switch plottype
    case 1
        for a = 1:length(labels)
            markerstyle = [colors(a) '.'];
            xpos1 = xpos(label == a);
            ypos1 = ypos(label == a);
            zpos1 = zpos(label == a);            
            plot3(xpos1,ypos1,zpos1,markerstyle);
            hold on;
        end
        hold off;
        axis('equal');
    case 2
        xpos1 = xpos(clusteridx>0);
        ypos1 = ypos(clusteridx>0);
        zpos1 = zpos(clusteridx>0);
        label1 = label(clusteridx>0);
        for a = 1:length(labels)
            markerstyle = [colors(a) '.'];
            xpos2 = xpos1(label1 == a);
            ypos2 = ypos1(label1 == a);
            zpos2 = zpos1(label1 == a);            
            plot3(xpos2,ypos2,zpos2,markerstyle);
            hold on;
        end
        hold off;
        axis('equal');
    case 3
        xpos1 = xpos(clusteridx==0);
        ypos1 = ypos(clusteridx==0);
        zpos1 = zpos(clusteridx==0);
        label1 = label(clusteridx==0);
        for a = 1:length(labels)
            markerstyle = [colors(a) '.'];
            xpos2 = xpos1(label1 == a);
            ypos2 = ypos1(label1 == a);
            zpos2 = zpos1(label1 == a);            
            plot3(xpos2,ypos2,zpos2,markerstyle);
            hold on;
        end
        hold off;
        axis('equal');
end

% --- Executes on selection change in paramselect.
function paramselect_Callback(hObject, eventdata, handles)
% hObject    handle to paramselect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns paramselect contents as cell array
%        contents{get(hObject,'Value')} returns selected item from paramselect


% --- Executes during object creation, after setting all properties.
function paramselect_CreateFcn(hObject, eventdata, handles)
% hObject    handle to paramselect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider2_Callback(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider



% --- Executes on button press in filtercluster.
function filtercluster_Callback(hObject, eventdata, handles)
% hObject    handle to filtercluster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isappdata(handles.clustertable1,'props1')
    props1 = getappdata(handles.clustertable1,'props1');
    shp1 = getappdata(handles.clustertable1,'shp1');
    datatab1 = handles.clustertable1.Data;
    prop1keep = min(props1 >= datatab1(:,1)' & props1 <= datatab1(:,2)',[],2);
    props1filt = props1(prop1keep,:);
    shp1filt = shp1(prop1keep);
    setappdata(handles.filtercluster,'props1filt',props1filt);
    setappdata(handles.filtercluster,'shp1filt',shp1filt);
    datatab1new = [min(props1filt)' max(props1filt)'];
    handles.clustertable1.Data = datatab1new;
    numclusters = length(shp1filt);
    axes(handles.axes1);
    for a = 1:numclusters
        plot(shp1filt{a},'FaceColor','r','FaceAlpha',0.5);
        hold on;
    end
end

if isappdata(handles.clustertable2,'props2')
    props2 = getappdata(handles.clustertable2,'props2');
    shp2 = getappdata(handles.clustertable2,'shp2');
    datatab2 = handles.clustertable2.Data;
    prop2keep = min(props2 >= datatab2(:,1)' & props2 <= datatab2(:,2)',[],2);
    props2filt = props2(prop2keep,:);
    shp2filt = shp2(prop2keep);
    setappdata(handles.filtercluster,'props2filt',props2filt);
    setappdata(handles.filtercluster,'shp2filt',shp2filt);
    datatab2new = [min(props2filt)' max(props2filt)'];
    handles.clustertable2.Data = datatab2new;
    numclusters = length(shp2filt);
    axes(handles.axes1);
    for a = 1:numclusters
        plot(shp2filt{a},'FaceColor','g','FaceAlpha',0.5);
        hold on;
    end
end

if isappdata(handles.clustertable3,'props3')
    props3 = getappdata(handles.clustertable3,'props3');
    shp3 = getappdata(handles.clustertable2,'shp3');
    datatab3 = handles.clustertable3.Data;
    prop3keep = min(props3 >= datatab3(:,1)' & props3 <= datatab3(:,2)',[],2);
    props3filt = props3(prop3keep,:);
    shp3filt = shp3(prop3keep);
    setappdata(handles.filtercluster,'props2filt',props3filt);
    setappdata(handles.filtercluster,'shp2filt',shp3filt);
    datatab3new = [min(props3filt)' max(props3filt)'];
    handles.clustertable3.Data = datatab3new;
    numclusters = length(shp3filt);
    axes(handles.axes1);
    for a = 1:numclusters
        plot(shp3filt{a},'FaceColor','b','FaceAlpha',0.5);
        hold on;
    end
end
hold off;


% --- Executes on button press in plothist.
function plothist_Callback(hObject, eventdata, handles)
% hObject    handle to plothist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
paramnum = handles.paramselect.Value;
plotcolor = ['r','g','b'];
logx = handles.logx.Value;
logy = handles.logy.Value;

if isappdata(handles.filtercluster,'props1filt')
    props1filt = getappdata(handles.filtercluster,'props1filt');
    [n1,x1] = histcounts(props1filt(:,paramnum));
    xc1 = mean([x1(1:(end-1)); x1(2:end)]);
    axes(handles.axes2);
    bar(xc1,n1,plotcolor(1));
    if logx
        handles.axes2.XScale = 'log';
    else
        handles.axes.XScale = 'linear';
    end
    if logy
        handles.axes2.YScale = 'log';
    else
        handles.axes2.YScale = 'linear';
    end
    hold on;
end

if isappdata(handles.filtercluster,'props2filt')
    props2filt = getappdata(handles.filtercluster,'props2filt');
    [n2,x2] = histcounts(props2filt(:,paramnum));
    xc2 = mean([x2(1:(end-1)); x2(2:end)]);
    axes(handles.axes2);
    bar(xc2,n2,plotcolor(2));
    if logx
        handles.axes2.XScale = 'log';
    else
        handles.axes.XScale = 'linear';
    end
    if logy
        handles.axes2.YScale = 'log';
    else
        handles.axes2.YScale = 'linear';
    end
    hold on;
end

if isappdata(handles.filtercluster,'props3filt')
    props3filt = getappdata(handles.filtercluster,'props3filt');
    [n3,x3] = histcounts(props3filt(:,paramnum));
    xc3 = mean([x3(1:(end-1)); x3(2:end)]);
    axes(handles.axes2);
    bar(xc3,n3,plotcolor(3));
    if logx
        handles.axes2.XScale = 'log';
    else
        handles.axes.XScale = 'linear';
    end
    if logy
        handles.axes2.YScale = 'log';
    else
        handles.axes2.YScale = 'linear';
    end
    hold on;
end
hold off;





% --- Executes on selection change in plottype.
function plottype_Callback(hObject, eventdata, handles)
% hObject    handle to plottype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns plottype contents as cell array
%        contents{get(hObject,'Value')} returns selected item from plottype


% --- Executes during object creation, after setting all properties.
function plottype_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plottype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in clustercalc.
function clustercalc_Callback(hObject, eventdata, handles)
% hObject    handle to clustercalc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data = getappdata(handles.loaddata,'data');
xpos = data(:,1);
ypos = data(:,2);
zpos = data(:,3);
label = data(:,4);
labels = unique(label);
clusteridx = data(:,5);
%clusterpos = cell(3,1);
axes(handles.axes1); hold on;
if ismember(1,labels)
    clusterpos1 = [xpos(label == 1 & clusteridx > 0) ypos(label == 1 & clusteridx > 0) zpos(label == 1 & clusteridx > 0) clusteridx(label == 1 & clusteridx > 0)];
    currclusters = unique(clusteridx(label == 1 & clusteridx > 0));
    numclusters = length(currclusters);
    if numclusters>0
        shp1 = cell(numclusters,1);
        %props1 = zeros(numclusters,1);
        for a = 1:numclusters
            currclusteridx = currclusters(a);
            currpos = clusterpos1(:,4) == currclusteridx;
            xyz = clusterpos1(currpos,1:3);
            x = xyz(:,1); y = xyz(:,2); z = xyz(:,3);
            shp11 = alphaShape(x,y,z);
            alpha = shp11.Alpha;
            shp1{a} = alphaShape(x,y,z,1.25*alpha);
            plot(shp1{a},'FaceColor','r','FaceAlpha',0.5);
            props1(a).numpoints = length(x);
            props1(a).centx = mean(x);
            props1(a).centy = mean(y);
            props1(a).centz = mean(z);
            props1(a).volume = volume(shp1{a});
            props1(a).SA = surfaceArea(shp1{a});
            props1(a).SAVR = surfaceArea(shp1{a})/volume(shp1{a});
            props1(a).deff = 6*volume(shp1{a})/surfaceArea(shp1{a});
        end
        %setappdata(handles.clustertable1,'props1',props1);
        setappdata(handles.clustertable1,'shp1',shp1);
        setappdata(handles.clustertable1,'params',fieldnames(props1));
        handles.paramselect.String = fieldnames(props1);
        props1mat = squeeze(cell2mat(struct2cell(props1)))';
        setappdata(handles.clustertable1,'props1',props1mat);
        mins1 = min(props1mat,[],1);
        maxs1 = max(props1mat,[],1);
        datatab1 = [mins1' maxs1'];
        handles.clustertable1.Data = datatab1;
        handles.clustertable1.RowName = fieldnames(props1);
        handles.clustertable1.ColumnName = {'Min' 'Max'};
        handles.clustertable1.ColumnEditable = logical([1 1]);
        setappdata(handles.filtercluster,'props1filt',props1mat);
        setappdata(handles.filtercluster,'shp1filt',shp1);
    end
end
        

if ismember(2,labels)
    clusterpos2 = [xpos(label == 2 & clusteridx > 0) ypos(label == 2 & clusteridx > 0) zpos(label == 2 & clusteridx > 0) clusteridx(label == 2 & clusteridx > 0)];
    currclusters = unique(clusteridx(label == 2 & clusteridx > 0));
    numclusters = length(currclusters);
    if numclusters>0
        shp2 = cell(numclusters,1);
        %props1 = zeros(numclusters,1);
        for a = 1:numclusters
            currclusteridx = currclusters(a);
            currpos = clusterpos2(:,4) == currclusteridx;
            xyz = clusterpos2(currpos,1:3);
            x = xyz(:,1); y = xyz(:,2); z = xyz(:,3);
            shp22 = alphaShape(x,y,z);
            alpha = shp22.Alpha;
            shp2{a} = alphaShape(x,y,z,1.25*alpha);
            plot(shp2{a},'FaceColor','g','FaceAlpha',0.5);
            props2(a).numpoints = length(x);
            props2(a).centx = mean(x);
            props2(a).centy = mean(y);
            props2(a).centz = mean(z);
            props2(a).volume = volume(shp2{a});
            props2(a).SA = surfaceArea(shp2{a});
            props2(a).SAVR = surfaceArea(shp2{a})/volume(shp2{a});
            props2(a).deff = 6*volume(shp2{a})/surfaceArea(shp2{a});
        end
        %setappdata(handles.clustertable2,'props2',props2);
        setappdata(handles.clustertable2,'shp2',shp2);
        setappdata(handles.clustertable2,'params',fieldnames(props2));
        props2mat = squeeze(cell2mat(struct2cell(props2)))';
        setappdata(handles.clustertable2,'props2',props2mat);
        mins2 = min(props2mat,[],1);
        maxs2 = max(props2mat,[],1);
        datatab2 = [mins2' maxs2'];
        handles.clustertable2.Data = datatab2;
        handles.clustertable2.RowName = fieldnames(props2);
        handles.clustertable2.ColumnName = {'Min' 'Max'};
        handles.clustertable2.ColumnEditable = logical([1 1]);
        setappdata(handles.filtercluster,'props2filt',props2mat);
        setappdata(handles.filtercluster,'shp2filt',shp2);
    end
end

if ismember(3,labels)
    clusterpos3 = [xpos(label == 3 & clusteridx > 0) ypos(label == 3 & clusteridx > 0) zpos(label == 3 & clusteridx > 0) clusteridx(label == 3 & clusteridx > 0)];
    currclusters = unique(clusteridx(label == 3 & clusteridx > 0));
    numclusters = length(currclusters);
    if numclusters>0
        shp3 = cell(numclusters,1);
        %props1 = zeros(numclusters,1);
        for a = 1:numclusters
            currclusteridx = currclusters(a);
            currpos = clusterpos3(:,4) == currclusteridx;
            xyz = clusterpos3(currpos,1:3);
            x = xyz(:,1); y = xyz(:,2); z = xyz(:,3);
            shp33 = alphaShape(x,y,z);
            alpha = shp33.Alpha;
            shp3{a} = alphaShape(x,y,z,1.25*alpha);
            plot(shp3{a},'FaceColor','b','FaceAlpha',0.5);
            props3(a).numpoints = length(x);
            props3(a).centx = mean(x);
            props3(a).centy = mean(y);
            props3(a).centz = mean(z);
            props3(a).volume = volume(shp3{a});
            props3(a).SA = surfaceArea(shp3{a});
            props3(a).SAVR = surfaceArea(shp3{a})/volume(shp3{a});
            props3(a).deff = 6*volume(shp3{a})/surfaceArea(shp3{a});
        end
        %setappdata(handles.clustertable3,'props3',props3);
        setappdata(handles.clustertable3,'shp3',shp3);
        setappdata(handles.clustertable3,'params',fieldnames(props3));
        props3mat = squeeze(cell2mat(struct2cell(props3)))';
        setappdata(handles.clustertable3,'props3',props3mat);
        mins3 = min(props3mat,[],1);
        maxs3 = max(props3mat,[],1);
        datatab3 = [mins3' maxs3'];
        handles.clustertable3.Data = datatab3;
        handles.clustertable3.RowName = fieldnames(props3);
        handles.clustertable3.ColumnName = {'Min' 'Max'};
        handles.clustertable3.ColumnEditable = logical([1 1]);
        setappdata(handles.filtercluster,'props3filt',props3mat);
        setappdata(handles.filtercluster,'shp3filt',shp3);
    end
end
hold off;



% --- Executes on selection change in distfrom.
function distfrom_Callback(hObject, eventdata, handles)
% hObject    handle to distfrom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns distfrom contents as cell array
%        contents{get(hObject,'Value')} returns selected item from distfrom


% --- Executes during object creation, after setting all properties.
function distfrom_CreateFcn(hObject, eventdata, handles)
% hObject    handle to distfrom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in distto.
function distto_Callback(hObject, eventdata, handles)
% hObject    handle to distto (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns distto contents as cell array
%        contents{get(hObject,'Value')} returns selected item from distto


% --- Executes during object creation, after setting all properties.
function distto_CreateFcn(hObject, eventdata, handles)
% hObject    handle to distto (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in calcdist.
function calcdist_Callback(hObject, eventdata, handles)
% hObject    handle to calcdist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
distfromparam = handles.distfrom.Value;
disttoparam = handles.distto.Value;

if distfromparam == 1
    datapeaks = getappdata(handles.loaddata,'data');
    labels = datapeaks(:,4);
    datafrom = datapeaks(labels==1,1:3);
elseif distfromparam == 2
    dataclusters = getappdata(handles.filtercluster,'props1filt');
    params = handles.clustertable1.RowName;
    xposidx = find(contains(params,'centx'),1);
    yposidx = find(contains(params,'centy'),1);
    zposidx = find(contains(params,'centz'),1);
    datafrom = dataclusters(:,[xposidx yposidx zposidx]);
elseif distfromparam == 3
    datapeaks = getappdata(handles.loaddata,'data');
    labels = datapeaks(:,4);
    datafrom = datapeaks(labels==2,1:3);
elseif distfromparam == 4
    dataclusters = getappdata(handles.filtercluster,'props2filt');
    params = handles.clustertable2.RowName;
    xposidx = find(contains(params,'centx'),1);
    yposidx = find(contains(params,'centy'),1);
    zposidx = find(contains(params,'centz'),1);
    datafrom = dataclusters(:,[xposidx yposidx zposidx]);
elseif distfromparam == 5
    datapeaks = getappdata(handles.loaddata,'data');
    labels = datapeaks(:,4);
    datafrom = datapeaks(labels==3,1:3);
elseif distfromparam == 6
    dataclusters = getappdata(handles.filtercluster,'props3filt');
    params = handles.clustertable3.RowName;
    xposidx = find(contains(params,'centx'),1);
    yposidx = find(contains(params,'centy'),1);
    zposidx = find(contains(params,'centz'),1);
    datafrom = dataclusters(:,[xposidx yposidx zposidx]);
end

if disttoparam == 1
    datapeaks = getappdata(handles.loaddata,'data');
    labels = datapeaks(:,4);
    datato = datapeaks(labels==1,1:3);
elseif disttoparam == 2
    dataclusters = getappdata(handles.filtercluster,'props1filt');
    params = handles.clustertable1.RowName;
    xposidx = find(contains(params,'centx'),1);
    yposidx = find(contains(params,'centy'),1);
    zposidx = find(contains(params,'centz'),1);
    datato = dataclusters(:,[xposidx yposidx zposidx]);
elseif disttoparam == 3
    datapeaks = getappdata(handles.loaddata,'data');
    labels = datapeaks(:,4);
    datato = datapeaks(labels==2,1:3);
elseif disttoparam == 4
    dataclusters = getappdata(handles.filtercluster,'props2filt');
    params = handles.clustertable2.RowName;
    xposidx = find(contains(params,'centx'),1);
    yposidx = find(contains(params,'centy'),1);
    zposidx = find(contains(params,'centz'),1);
    datato = dataclusters(:,[xposidx yposidx zposidx]);
elseif disttoparam == 5
    datapeaks = getappdata(handles.loaddata,'data');
    labels = datapeaks(:,4);
    datato = datapeaks(labels==3,1:3);
elseif disttoparam == 6
    dataclusters = getappdata(handles.filtercluster,'props3filt');
    params = handles.clustertable3.RowName;
    xposidx = find(contains(params,'centx'),1);
    yposidx = find(contains(params,'centy'),1);
    zposidx = find(contains(params,'centz'),1);
    datato = dataclusters(:,[xposidx yposidx zposidx]);
end

if distfromparam == disttoparam
    distmat = pdist(datafrom);
    distmat = squareform(distmat);
    distmat(distmat == 0) = NaN;
else
    distmat = pdist2(datafrom,datato,'euclidean');
end
nndist = min(distmat,[],2);
setappdata(handles.calcdist,'nndist',nndist);
axes(handles.axes3);
cla;
[n,x] = histcounts(nndist);
xc = mean([x(1:(end-1)); x(2:end)]);
bar(xc,n);
xlabel('Distance (nm)');
ylabel('Count');
mediannndist = median(nndist);
stdnndist = std(nndist);
txt = ['median: ' num2str(mediannndist) '+/-' num2str(stdnndist) ' nm'];
text(mean(x),mean(n),txt);

% --- Executes on button press in saveresults.
function saveresults_Callback(hObject, eventdata, handles)
% hObject    handle to saveresults (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
output.peaksdata = getappdata(handles.loaddata,'data');
try
    output.Label1.shapes = getappdata(handles.clustertable1,'shp1');
    output.Label1.clusterprops = getappdata(handles.clustertable1,'props1');
    output.Label1.parameternames = handles.clustertable1.RowName;
    output.Label1.paramfilters = handles.clustertable1.Data;
    output.Label1.filteredclusters = getappdata(handles.filtercluster,'props1filt');
    output.Label1.filteredshapes = getappdata(handles.filtercluster,'shp1filt');
catch
end

try
    output.Label2.shapes = getappdata(handles.clustertable2,'shp2');
    output.Label2.clusterprops = getappdata(handles.clustertable2,'props2');
    output.Label2.parameternames = handles.clustertable2.RowName;
    output.Label2.paramfilters = handles.clustertable2.Data;
    output.Label2.filteredclusters = getappdata(handles.filtercluster,'props2filt');
    output.Label2.filteredshapes = getappdata(handles.filtercluster,'shp2filt');
catch
end

try
    output.Label3.shapes = getappdata(handles.clustertable3,'shp3');
    output.Label3.clusterprops = getappdata(handles.clustertable3,'props3');
    output.Label3.parameternames = handles.clustertable3.RowName;
    output.Label3.paramfilters = handles.clustertable3.Data;
    output.Label3.filteredclusters = getappdata(handles.filtercluster,'props3filt');
    output.Label3.filteredshapes = getappdata(handles.filtercluster,'shp3filt');
catch
end

diststring = handles.distfrom.String;
distfromparam = handles.distfrom.Value;
disttoparam = handles.distto.Value;
output.NNdistFrom = diststring(distfromparam);
output.NNdistTo = diststring(disttoparam);
output.NNdistances = getappdata(handles.calcdist,'nndist');
fullpath = handles.currfile.String;
[path,filename] = fileparts(fullpath);
[filesave,pathsave] = uiputfile('*.mat','Save Full Results As...',fullfile(path,[filename '.mat']));
save(fullfile(pathsave,filesave),'output');




% --- Executes on button press in reload.
function reload_Callback(hObject, eventdata, handles)
% hObject    handle to reload (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isappdata(handles.clustertable1,'props1')
    props1 = getappdata(handles.clustertable1,'props1');
    shp1 = getappdata(handles.clustertable1,'shp1');
    datatab1 = [min(props1)' max(props1)'];
    handles.clustertable1.Data = datatab1;
    setappdata(handles.filtercluster,'props1filt',props1);
    setappdata(handles.filtercluster,'shp1filt',shp1);
end

if isappdata(handles.clustertable2,'props2')
    props2 = getappdata(handles.clustertable2,'props2');
    shp2 = getappdata(handles.clustertable2,'shp2');
    datatab2 = [min(props2)' max(props2)'];
    handles.clustertable2.Data = datatab2;
    setappdata(handles.filtercluster,'props2filt',props2);
    setappdata(handles.filtercluster,'shp2filt',shp2);
end

if isappdata(handles.clustertable3,'props3')
    props3 = getappdata(handles.clustertable3,'props3');
    shp3 = getappdata(handles.clustertable3,'shp3');
    datatab3 = [min(props3)' max(props3)'];
    handles.clustertable3.Data = datatab3;
    setappdata(handles.filtercluster,'props3filt',props3);
    setappdata(handles.filtercluster,'shp3filt',shp3);
end

% --- Executes on button press in logx.
function logx_Callback(hObject, eventdata, handles)
% hObject    handle to logx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of logx


% --- Executes on button press in logy.
function logy_Callback(hObject, eventdata, handles)
% hObject    handle to logy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of logy
