function varargout = ipalmplotter(varargin)
% IPALMPLOTTER MATLAB code for ipalmplotter.fig
%      IPALMPLOTTER, by itself, creates a new IPALMPLOTTER or raises the existing
%      singleton*.
%
%      H = IPALMPLOTTER returns the handle to a new IPALMPLOTTER or the handle to
%      the existing singleton*.
%
%      IPALMPLOTTER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IPALMPLOTTER.M with the given input arguments.
%
%      IPALMPLOTTER('Property','Value',...) creates a new IPALMPLOTTER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ipalmplotter_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ipalmplotter_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ipalmplotter

% Last Modified by GUIDE v2.5 11-Aug-2020 12:48:25

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ipalmplotter_OpeningFcn, ...
                   'gui_OutputFcn',  @ipalmplotter_OutputFcn, ...
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


% --- Executes just before ipalmplotter is made visible.
function ipalmplotter_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ipalmplotter (see VARARGIN)

% Choose default command line output for ipalmplotter
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ipalmplotter wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ipalmplotter_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in xvals.
function xvals_Callback(hObject, eventdata, handles)
% hObject    handle to xvals (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns xvals contents as cell array
%        contents{get(hObject,'Value')} returns selected item from xvals


% --- Executes during object creation, after setting all properties.
function xvals_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xvals (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in yvals.
function yvals_Callback(hObject, eventdata, handles)
% hObject    handle to yvals (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns yvals contents as cell array
%        contents{get(hObject,'Value')} returns selected item from yvals


% --- Executes during object creation, after setting all properties.
function yvals_CreateFcn(hObject, eventdata, handles)
% hObject    handle to yvals (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in zvals.
function zvals_Callback(hObject, eventdata, handles)
% hObject    handle to zvals (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns zvals contents as cell array
%        contents{get(hObject,'Value')} returns selected item from zvals


% --- Executes during object creation, after setting all properties.
function zvals_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zvals (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in histvals.
function histvals_Callback(hObject, eventdata, handles)
% hObject    handle to histvals (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns histvals contents as cell array
%        contents{get(hObject,'Value')} returns selected item from histvals


% --- Executes during object creation, after setting all properties.
function histvals_CreateFcn(hObject, eventdata, handles)
% hObject    handle to histvals (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
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


% --- Executes on button press in plothist.
function plothist_Callback(hObject, eventdata, handles)
% hObject    handle to plothist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.histtable.Data = [];
paramnum = get(handles.histvals,'Value');
datafilt = getappdata(handles.filter,'datafilt');
labels = unique(datafilt(:,27));
numlab = length(unique(datafilt(:,27)));
histdata = cell(numlab,1);
for a = 1:numlab
    labindx = datafilt(:,27) == labels(a);
    histdata{a,1} = datafilt(labindx,paramnum);
    if min(histdata{a,1})<=0
        set(handles.logx,'Value',0);
    end
end

logx = get(handles.logx,'Value');
logy = get(handles.logy,'Value');
n = cell(numlab,1);
x = cell(numlab,1);
plotcolor = {'r','g','b'};
peakval = cell(numlab,1);
meanval = cell(numlab,1);
medval = cell(numlab,1);
fwhmval = cell(numlab,1);
gaussfit = get(handles.gaussfit,'Value');
gaussamp = cell(numlab,1);
gausscenter = cell(numlab,1);
gausssigma = cell(numlab,1);
axes(handles.axes2);
for a = 1:numlab
    if logx
        [n{a},x{a}] = histcounts(log10(histdata{a,1}));
    else
        [n{a},x{a}] = histcounts(histdata{a,1});
    end
    
    if logy
        n{a} = log10(n{a});
    end
    cents = mean([x{a}(1:(end-1));x{a}(2:end)],1);
    bar(cents,n{a},plotcolor{a}); hold on
    [~,peakind] = max(n{a});
    peakval{a} = x{a}(peakind); handles.histtable.Data{1,a} = peakval{a};
    meanval{a} = mean(histdata{a}); handles.histtable.Data{2,a} = meanval{a};
    medval{a} = median(histdata{a}); handles.histtable.Data{3,a} = medval{a};
    try
        fwhmval{a} = fwhm(x{a},n{a}); handles.histtable.Data{4,a} = fwhmval{a};
    catch
        fwhmval{a} = 'Not found';
    end
    if gaussfit
        ft = fittype('gauss1');
        try
            fm = fit(cents',n{a}',ft);
            plot(fm,plotcolor{a})
            coeffs = coeffvalues(fm);
            gaussamp{a} = coeffs(1); 
            gausscenter{a} = coeffs(2);
            gausssigma{a} = coeffs(3)/2;            
        catch
            gaussamp{a} = 'Bad Fit';
            gausscenter{a} = 'Bad Fit';
            gausssigma{a} = 'Bad Fit'; 
        end
        handles.histtable.Data{5,a} = gaussamp{a};
        handles.histtable.Data{6,a} = gausscenter{a};
        handles.histtable.Data{7,a} = gausssigma{a};
    end
end
legend;
hold off
histvalues.bins = x;
histvalues.counts = n;
exporthistdata.histvaules = histvalues;
params = get(handles.histvals,'String');
summarydata.parametername = params(paramnum);
summarydata.peakvalue = handles.histtable.Data(1,:);
summarydata.meanvalue = handles.histtable.Data(2,:);
summarydata.medianvalue = handles.histtable.Data(3,:);
summarydata.fwhm = handles.histtable.Data(4,:);
if gaussfit
    summarydata.gaussamp = handles.histtable.Data(5,:);
    summarydata.gausscenter = handles.histtable.Data(6,:);
    summarydata.gausssigma = handles.histtable.Data(7,:);
else
    summarydata.gaussamp = NaN;
    summarydata.gausscenter = NaN;
    summarydata.gausssigma = NaN;
end
    
exporthistdata.summarydata = summarydata;
setappdata(handles.histtable,'exporthistdata',exporthistdata);


% --- Executes on button press in plotvals.
function plotvals_Callback(hObject, eventdata, handles)
% hObject    handle to plotvals (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

datafilt = getappdata(handles.filter,'datafilt');
nparts = size(datafilt,1);
handles.numparticles.String = num2str(nparts);

xparamval = get(handles.xvals,'Value');
xparamlist = get(handles.xvals,'String');
xparamname = xparamlist{xparamval};
xdata = datafilt(:,xparamval);

yparamval = get(handles.yvals,'Value');
yparamlist = get(handles.yvals,'String');
yparamname = yparamlist{yparamval};
ydata = datafilt(:,yparamval);

zparamval = get(handles.zvals,'Value');
if zparamval > 1
    zparamlist = get(handles.zvals,'String');
    zparamname = zparamlist{zparamval};
    zdata = datafilt(:,zparamval-1);
end

markerstyleval = handles.markerstyle.Value;
markerstylelist = handles.markerstyle.String;
markerstyle = markerstylelist{markerstyleval};
markersize = str2double(handles.markersize.String);

markercolor = handles.spotcolor.Value;

if zparamval == 1 && markercolor == 1
    plottype = 1; %2d plot no color
elseif zparamval == 1 && markercolor == 2
    plottype = 2; %2d plot, color by label 
elseif zparamval == 1 && markercolor == 3
    plottype = 1; %2d plot, color by z
elseif zparamval > 1 && markercolor == 1
    plottype = 3; %3d plot, no color
elseif zparamval > 1 && markercolor == 2
    plottype = 4; %3d plot color by label
elseif zparamval > 1 && markercolor == 3
    plottype = 5; %3d plot color by z
end

axes(handles.axes1);
cla

switch plottype
    case 1
        marker = [markerstyle 'k'];
        plot(xdata,ydata,marker,'MarkerSize',markersize);
        xlabel(xparamname);
        ylabel(yparamname);
        hold off;
    case 2
        markercolor = ['r','g','b'];
        labels = unique(datafilt(:,27));
        numlabels = length(labels);
        xdata2 = cell(numlabels,1);
        ydata2 = cell(numlabels,1);
        marker = cell(numlabels,1);
        for a = 1:numlabels
            labelindex = datafilt(:,27) == labels(a);
            xdata2{a} = xdata(labelindex); 
            ydata2{a} = ydata(labelindex); 
            marker{a} = [markerstyle markercolor(label(a))]; 
            plot(xdata2{a},ydata2{a},marker{a},'MarkerSize',markersize);
            hold on;
        end
            xlabel(xparamname);
            ylabel(yparamname);
            hold off
    case 3
        marker = [markerstyle 'k'];
        plot3(xdata,ydata,zdata,marker,'MarkerSize',markersize);
        if handles.axisequal.Value
            axis('equal');
        end
        xlabel(xparamname);
        ylabel(yparamname);
        zlabel(zparamname);
        hold off
    
    case 4
        markercolor = ['r','g','b'];
        labels = unique(datafilt(:,27));
        numlabels = length(labels);
        xdata2 = cell(numlabels,1);
        ydata2 = cell(numlabels,1);
        zdata2 = cell(numlabels,1);
        marker = cell(numlabels,1);
        for a = 1:numlabels
            labelindex = datafilt(:,27) == labels(a);
            xdata2{a} = xdata(labelindex); 
            ydata2{a} = ydata(labelindex);
            zdata2{a} = zdata(labelindex);
            marker{a} = [markerstyle markercolor(labels(a))]; 
            plot3(xdata2{a},ydata2{a},zdata2{a},marker{a},'MarkerSize',markersize);
            hold on;
        end
        if handles.axisequal.Value
            axis('equal');
        end
        xlabel(xparamname);
        ylabel(yparamname);
        zlabel(zparamname);
        hold off
    case 5
        map = [0.242200000000000,0.150400000000000,0.660300000000000;0.250390476190476,0.164995238095238,0.707614285714286;0.257771428571429,0.181780952380952,0.751138095238095;0.264728571428571,0.197757142857143,0.795214285714286;0.270647619047619,0.214676190476190,0.836371428571429;0.275114285714286,0.234238095238095,0.870985714285714;0.278300000000000,0.255871428571429,0.899071428571429;0.280333333333333,0.278233333333333,0.922100000000000;0.281338095238095,0.300595238095238,0.941376190476191;0.281014285714286,0.322757142857143,0.957885714285714;0.279466666666667,0.344671428571429,0.971676190476191;0.275971428571429,0.366680952380952,0.982904761904762;0.269914285714286,0.389200000000000,0.990600000000000;0.260242857142857,0.412328571428571,0.995157142857143;0.244033333333333,0.435833333333333,0.998833333333333;0.220642857142857,0.460257142857143,0.997285714285714;0.196333333333333,0.484719047619048,0.989152380952381;0.183404761904762,0.507371428571429,0.979795238095238;0.178642857142857,0.528857142857143,0.968157142857143;0.176438095238095,0.549904761904762,0.952019047619048;0.168742857142857,0.570261904761905,0.935871428571429;0.154000000000000,0.590200000000000,0.921800000000000;0.146028571428571,0.609119047619048,0.907857142857143;0.138023809523810,0.627628571428572,0.897290476190476;0.124814285714286,0.645928571428571,0.888342857142857;0.111252380952381,0.663500000000000,0.876314285714286;0.0952095238095238,0.679828571428571,0.859780952380952;0.0688714285714285,0.694771428571429,0.839357142857143;0.0296666666666667,0.708166666666667,0.816333333333333;0.00357142857142858,0.720266666666667,0.791700000000000;0.00665714285714287,0.731214285714286,0.766014285714286;0.0433285714285715,0.741095238095238,0.739409523809524;0.0963952380952380,0.750000000000000,0.712038095238095;0.140771428571429,0.758400000000000,0.684157142857143;0.171700000000000,0.766961904761905,0.655442857142857;0.193766666666667,0.775766666666667,0.625100000000000;0.216085714285714,0.784300000000000,0.592300000000000;0.246957142857143,0.791795238095238,0.556742857142857;0.290614285714286,0.797290476190476,0.518828571428572;0.340642857142857,0.800800000000000,0.478857142857143;0.390900000000000,0.802871428571429,0.435447619047619;0.445628571428572,0.802419047619048,0.390919047619048;0.504400000000000,0.799300000000000,0.348000000000000;0.561561904761905,0.794233333333333,0.304480952380953;0.617395238095238,0.787619047619048,0.261238095238095;0.671985714285714,0.779271428571429,0.222700000000000;0.724200000000000,0.769842857142857,0.191028571428571;0.773833333333333,0.759804761904762,0.164609523809524;0.820314285714286,0.749814285714286,0.153528571428571;0.863433333333333,0.740600000000000,0.159633333333333;0.903542857142857,0.733028571428571,0.177414285714286;0.939257142857143,0.728785714285714,0.209957142857143;0.972757142857143,0.729771428571429,0.239442857142857;0.995647619047619,0.743371428571429,0.237147619047619;0.996985714285714,0.765857142857143,0.219942857142857;0.995204761904762,0.789252380952381,0.202761904761905;0.989200000000000,0.813566666666667,0.188533333333333;0.978628571428571,0.838628571428572,0.176557142857143;0.967647619047619,0.863900000000000,0.164290476190476;0.961009523809524,0.889019047619048,0.153676190476191;0.959671428571429,0.913457142857143,0.142257142857143;0.962795238095238,0.937338095238095,0.126509523809524;0.969114285714286,0.960628571428571,0.106361904761905;0.976900000000000,0.983900000000000,0.0805000000000000];
        zslices = linspace(min(zdata),max(zdata),length(map));
        for a = 1:(length(zslices)-1)
            currz = zdata(zdata>=zslices(a) & zdata<zslices(a+1));
            currx = xdata(zdata>=zslices(a) & zdata<zslices(a+1));
            curry = ydata(zdata>=zslices(a) & zdata<zslices(a+1));
            plot3(currx,curry,currz,markerstyle,'Color',map(a,:));
            hold on
        end
        if handles.axisequal.Value
            axis('equal');
        end
        xlabel(xparamname);
        ylabel(yparamname);
        zlabel(zparamname);
        hold off
    otherwise
        disp('not yet');
end



% --- Executes on button press in loaddata.
function loaddata_Callback(hObject, eventdata, handles)
% hObject    handle to loaddata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[datamat,paramlist,filename,pathname] = ipalmimport(133);
for a = 1:length(paramlist)
    paramlist{a} = regexprep(paramlist{a},' ','_');
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%If loading unanalyzed data, add new parameter values for new point analysis
if length(paramlist)<=49
    numpoints = size(datamat,1);
    datamat(:,length(paramlist)+1) = zeros(numpoints,1);
    paramlist{end+1} = 'Cluster_Index';
    datamat(:,length(paramlist)+1) = zeros(numpoints,1);
    paramlist{end+1} = 'Coloc_Value';
    datamat(:,length(paramlist)+1) = zeros(numpoints,1);
    paramlist{end+1} = 'NN_Distance';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%setappdata(handles.loaddata,'datafull',datafull);
setappdata(handles.loaddata,'paramlist',paramlist);
set(handles.filepath,'String',fullfile(pathname,filename));
set(handles.xvals,'String',paramlist);
set(handles.yvals,'String',paramlist);
set(handles.histvals,'String',paramlist);
paramlistz = cell(length(paramlist)+1,1);
paramlistz(1) = {'None'};
paramlistz(2:end) = paramlist;
set(handles.zvals,'String',paramlistz);
labelidx = find(contains(paramlist,'Label_Set'),1);
uniquelabs = unique(datamat(:,labelidx));
if uniquelabs == 0
    datamat(:,labelidx) = 1;
end
mins = min(datamat)';
maxs = max(datamat)';
datatab = [mins maxs];
handles.datatable.Data = datatab;
handles.datatable.RowName = paramlist;
handles.datatable.ColumnName = {'Min' 'Max'};
handles.datatable.ColumnEditable = logical([1 1]);
setappdata(handles.filter,'datamat',datamat);
setappdata(handles.filter,'datafilt',datamat);
setappdata(handles.unzoom,'datatab',datatab);
numpts = size(datamat,1);
handles.numparticles.String = num2str(numpts);

minindx = datamat >= datatab(:,1)';
maxindx = datamat <= datatab(:,2)';
keepindx = minindx & maxindx;
keepindx2 = min(keepindx,[],2)';
setappdata(handles.filter,'keepindx',keepindx2);


% --- Executes on selection change in spotcolor.
function spotcolor_Callback(hObject, eventdata, handles)
% hObject    handle to spotcolor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns spotcolor contents as cell array
%        contents{get(hObject,'Value')} returns selected item from spotcolor


% --- Executes during object creation, after setting all properties.
function spotcolor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to spotcolor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in filter.
function filter_Callback(hObject, eventdata, handles)
% hObject    handle to filter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
datamat = getappdata(handles.filter,'datamat');
datatab = handles.datatable.Data;
minindx = datamat >= datatab(:,1)';
maxindx = datamat <= datatab(:,2)';
keepindx = minindx & maxindx;
keepindx2 = min(keepindx,[],2)';
datafilt = datamat(keepindx2,:);
setappdata(handles.filter,'datafilt',datafilt);
setappdata(handles.unzoom,'datatab',datatab);
numpts = size(datafilt,1);
handles.numparticles.String = num2str(numpts);
setappdata(handles.filter,'keepindx',keepindx2);


% --- Executes on button press in reload.
function reload_Callback(hObject, eventdata, handles)
% hObject    handle to reload (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
datamat = getappdata(handles.filter,'datamat');
mins = min(datamat)';
maxs = max(datamat)';
datatab = [mins maxs];
handles.datatable.Data = datatab;
setappdata(handles.filter,'datafilt',datamat);
setappdata(handles.unzoom,'datatab',datatab);
numpts = size(datamat,1);
handles.numparticles.String = num2str(numpts);
minindx = datamat >= datatab(:,1)';
maxindx = datamat <= datatab(:,2)';
keepindx = minindx & maxindx;
keepindx2 = min(keepindx,[],2)';
setappdata(handles.filter,'keepindx',keepindx2);


% --- Executes on button press in gaussfit.
function gaussfit_Callback(hObject, eventdata, handles)
% hObject    handle to gaussfit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of gaussfit


% --- Executes on button press in exporthist.
function exporthist_Callback(hObject, eventdata, handles)
% hObject    handle to exporthist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
exporthistdata = getappdata(handles.histtable,'exporthistdata');
[filesave,pathsave] = uiputfile('*.mat','Save Data as (.mat)');
save(fullfile(pathsave,filesave),'exporthistdata');
 


% --- Executes on selection change in markerstyle.
function markerstyle_Callback(hObject, eventdata, handles)
% hObject    handle to markerstyle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns markerstyle contents as cell array
%        contents{get(hObject,'Value')} returns selected item from markerstyle


% --- Executes during object creation, after setting all properties.
function markerstyle_CreateFcn(hObject, eventdata, handles)
% hObject    handle to markerstyle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function markersize_Callback(hObject, eventdata, handles)
% hObject    handle to markersize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of markersize as text
%        str2double(get(hObject,'String')) returns contents of markersize as a double


% --- Executes during object creation, after setting all properties.
function markersize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to markersize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in axisequal.
function axisequal_Callback(hObject, eventdata, handles)
% hObject    handle to axisequal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of axisequal


% --- Executes on button press in zoomrect.
function zoomrect_Callback(hObject, eventdata, handles)
% hObject    handle to zoomrect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
datatab = handles.datatable.Data; 
axes(handles.axes1);
revz = handles.revz.Value;
if ~revz
    view([0 90])
elseif revz
    view([0 -90])
end
xvalnum = handles.xvals.Value;
yvalnum = handles.yvals.Value;
rect = getrect;
xvalmin = rect(1); 
yvalmin = rect(2);
xvalmax = rect(1)+rect(3);
yvalmax = rect(2)+rect(4);
datatab(xvalnum,1) = xvalmin;
datatab(xvalnum,2) = xvalmax;
datatab(yvalnum,1) = yvalmin;
datatab(yvalnum,2) = yvalmax;
handles.datatable.Data = datatab;
datamat = getappdata(handles.filter,'datamat');
minindx = datamat >= datatab(:,1)';
maxindx = datamat <= datatab(:,2)';
keepindx = minindx & maxindx;
keepindx2 = min(keepindx,[],2)';
datafilt = datamat(keepindx2,:);
setappdata(handles.filter,'datafilt',datafilt);
plotvals_Callback(handles.plotvals, eventdata, handles);
axes(handles.axes1);
if ~revz
    view([0 90])
elseif revz
    view([0 -90])
end
setappdata(handles.filter,'keepindx',keepindx2);


% --- Executes on button press in zoompoly.
function zoompoly_Callback(hObject, eventdata, handles)
% hObject    handle to zoompoly (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
datatab = handles.datatable.Data;
datamat = getappdata(handles.filter,'datamat');
axes(handles.axes1);
revz = handles.revz.Value;
if ~revz
    view([0 90])
elseif revz
    view([0 -90])
end
xvalnum = handles.xvals.Value;
yvalnum = handles.yvals.Value;
[xi,yi] = getline(handles.axes1,'closed');
xq = datamat(:,xvalnum);
yq = datamat(:,yvalnum);
in = inpolygon(xq,yq,xi,yi);
datapolyfilt = datamat(in,:);
minindx = datapolyfilt >= datatab(:,1)';
maxindx = datapolyfilt <= datatab(:,2)';
keepindx = minindx & maxindx;
keepindx2 = min(keepindx,[],2)';
datafilt = datapolyfilt(keepindx2,:);
setappdata(handles.filter,'keepindx',keepindx2);
mins = min(datafilt)';
maxs = max(datafilt)';
handles.datatable.Data = [mins maxs];
setappdata(handles.filter,'datafilt',datafilt);
plotvals_Callback(handles.plotvals, eventdata, handles);
axes(handles.axes1);
if ~revz
    view([0 90])
elseif revz
    view([0 -90])
end


% --- Executes on button press in unzoom.
function unzoom_Callback(hObject, eventdata, handles)
% hObject    handle to unzoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
datatab = getappdata(handles.unzoom,'datatab');
datamat = getappdata(handles.filter,'datamat');
set(handles.datatable,'Data',datatab);
minindx = datamat >= datatab(:,1)';
maxindx = datamat <= datatab(:,2)';
keepindx = minindx & maxindx;
keepindx2 = min(keepindx,[],2)';
datafilt = datamat(keepindx2,:);
setappdata(handles.filter,'keepindx',keepindx2);
setappdata(handles.filter,'datafilt',datafilt);
plotvals_Callback(handles.plotvals, eventdata, handles);
axes(handles.axes1);
revz = handles.revz.Value;
if ~revz
    view([0 90]);
elseif revz
    view([0 -90]);
end

% --- Executes on button press in exportdata.
function exportdata_Callback(hObject, eventdata, handles)
% hObject    handle to exportdata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
datafilt = getappdata(handles.filter,'datafilt');
paramlist = getappdata(handles.loaddata,'paramlist');
output = inputdlg('Pixel Size','Specify pixel size in nm.  Type 1 to save x/y pos in nm units',1,{'133'});
pixelsize = str2double(cell2mat(output));
xposidx = find(contains(paramlist,'X_Position'),1);
yposidx = find(contains(paramlist,'Y_Position'),1);
gxposidx = find(contains(paramlist,'Group X_Position'),1);
gyposidx = find(contains(paramlist,'Group Y_Position'),1);
datafilt(:,xposidx) = datafilt(:,xposidx)/pixelsize;
datafilt(:,yposidx) = datafilt(:,yposidx)/pixelsize;
datafilt(:,gxposidx) = datafilt(:,gxposidx)/pixelsize;
datafilt(:,gyposidx) = datafilt(:,gyposidx)/pixelsize;
[filesave,pathsave] = uiputfile('*.txt','Save Data as ASCII text...');
fullsave = fullfile(pathsave,filesave);
fid = fopen(fullsave,'wt');
for a = 1:length(paramlist)
    if a<length(paramlist)
        fprintf(fid,'%s\t',paramlist{a});
    elseif a==length(paramlist)
        fprintf(fid,'%s\n',paramlist{a});
    end
end
fclose(fid);       
dlmwrite(fullsave,datafilt,'delimiter','\t','-append');


% --- Executes on button press in xzproj.
function xzproj_Callback(hObject, eventdata, handles)
% hObject    handle to xzproj (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes1);
view([0 0])


% --- Executes on button press in yzproj.
function yzproj_Callback(hObject, eventdata, handles)
% hObject    handle to yzproj (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes1);
view([90 0])


% --- Executes on button press in xyproj.
function xyproj_Callback(hObject, eventdata, handles)
% hObject    handle to xyproj (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
revz = handles.revz.Value;
axes(handles.axes1);
if ~revz
    view([0 90])
elseif revz
    view([0 -90])
end
    


% --- Executes when entered data in editable cell(s) in datatable.
function datatable_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to datatable (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in hist2d.
function hist2d_Callback(hObject, eventdata, handles)
% hObject    handle to hist2d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
datafilt = getappdata(handles.filter,'datafilt');
nbinsx = str2double(handles.nbinsx.String);
nbinsy = str2double(handles.nbinsy.String);
nbins = [nbinsx nbinsy];
xparamnum = handles.xvals.Value;
yparamnum = handles.yvals.Value;
xparamnames = handles.xvals.String;
yparamnames = handles.yvals.String;
xname = xparamnames{xparamnum};
yname = yparamnames{yparamnum};
data = [datafilt(:,xparamnum) datafilt(:,yparamnum)];
axes(handles.axes1);
hist3(data,'Nbins',nbins,'CDataMode','auto','FaceColor','interp');
view([0 90])
axis('image');
xlabel(xname);
ylabel(yname);

function nbinsx_Callback(hObject, eventdata, handles)
% hObject    handle to nbinsx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nbinsx as text
%        str2double(get(hObject,'String')) returns contents of nbinsx as a double


% --- Executes during object creation, after setting all properties.
function nbinsx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nbinsx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nbinsy_Callback(hObject, eventdata, handles)
% hObject    handle to nbinsy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nbinsy as text
%        str2double(get(hObject,'String')) returns contents of nbinsy as a double


% --- Executes during object creation, after setting all properties.
function nbinsy_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nbinsy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in polysurf.
function polysurf_Callback(hObject, eventdata, handles)
% hObject    handle to polysurf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
datafilt = getappdata(handles.filter,'datafilt');
prompt = {'Polynomial Degree in X (1-5)','Polynomial Degree in Y (1-5)'};
title = 'Choose Polynomial Degree.  Use vector for multi-label.';
definput = {'1', '1'};
outs = inputdlg(prompt,title,definput,'on');




% --- Executes on button press in rundbscan.
function rundbscan_Callback(hObject, eventdata, handles)
% hObject    handle to rundbscan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
datafilt = getappdata(handles.filter,'datafilt');
datamat = getappdata(handles.filter,'datamat');
keepindx = getappdata(handles.filter,'keepindx');
paramlist = getappdata(handles.loaddata,'paramlist');
prompt = {'eps', 'Min. Points','Color Channel (1=red, 2=green, 3=blue)'};
dlg_title = 'Input DBSCAN parameters';
num_lines = 1;
defaultans = {'100','5','1'};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
epsilon = str2double(answer{1});
minpts = str2double(answer{2});
label = str2double(answer{3});
xposidx = find(contains(paramlist,'X_Position'),1);
yposidx = find(contains(paramlist,'Y_Position'),1);
zposidx = find(contains(paramlist,'Unwrapped_Z'),1);
labelidx = find(contains(paramlist,'Label_Set'),1);
xpos = datafilt(:,xposidx);
ypos = datafilt(:,yposidx);
zpos = datafilt(:,zposidx);
labkeep = datafilt(:,labelidx) == label;
xpos = xpos(labkeep);
ypos = ypos(labkeep);
zpos = zpos(labkeep);
P = [xpos ypos zpos];
indx = DBSCAN(P,epsilon,minpts);
plotvals_Callback(hObject, eventdata, handles);
set(handles.figure1,'CurrentAxes',handles.axes1); hold on;
for a = 1:max(indx)
    currindx = indx == a;
    currP = P(currindx,:);
    plot3(currP(:,1),currP(:,2),currP(:,3),'o','MarkerSize',8,'Color',0.5*rand(1,3)+0.5);       
end
view(0,90)

%Remove any previous DBSCAN results for the current label
clustercol = find(contains(paramlist,'Cluster_Index'),1);
fulllabel = datamat(:,labelidx);
fullabelidx = fulllabel==label;
datamat(fullabelidx,clustercol) = -inf;

%add in cluster index values in datafilt for current label
datafilt(labkeep,clustercol) = indx;
setappdata(handles.filter,'datafilt',datafilt);

%update filter table
datatab = handles.datatable.Data;
datatab(clustercol,1) = min(datafilt(:,clustercol));
datatab(clustercol,2) = max(datafilt(:,clustercol));
handles.datatable.Data= datatab;

%update "full" data on the current points/label
fullindex = keepindx' & fullabelidx;
datamat(fullindex,clustercol) = indx;
setappdata(handles.filter,'datamat',datamat);


% --- Executes on button press in exportanalysis.
function exportanalysis_Callback(hObject, eventdata, handles)
% hObject    handle to exportanalysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in coloc1.
function coloc1_Callback(hObject, eventdata, handles)
% hObject    handle to coloc1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%get data and user inputs
datafilt = getappdata(handles.filter,'datafilt');
datamat = getappdata(handles.filter,'datamat');
keepindx = getappdata(handles.filter,'keepindx');
paramlist = getappdata(handles.loaddata,'paramlist');
xposidx = find(contains(paramlist,'X_Position'),1);
yposidx = find(contains(paramlist,'Y_Position'),1);
zposidx = find(contains(paramlist,'Unwrapped_Z'),1);
labelidx = find(contains(paramlist,'Label_Set'),1);
colocidx = find(contains(paramlist,'Coloc_Value'),1);

prompt = {'Direction (1 = red to green, 2 = green to red','Max search radius (nm)','Search radius increment (nm)'};
dlg_title = 'Input coloc parameters';
num_lines = 1;
defaultans = {'1','100','1'};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
direction = str2double(answer{1});
rmax = str2double(answer{2});
rstep = str2double(answer{3});
r = rstep:rstep:rmax;

%extract xyz coords and separate by color
xpos = datafilt(:,xposidx);
ypos = datafilt(:,yposidx);
zpos = datafilt(:,zposidx);
xposr = xpos(datafilt(:,labelidx) == 1);
yposr = ypos(datafilt(:,labelidx) == 1);
zposr = zpos(datafilt(:,labelidx) == 1);
xposg = xpos(datafilt(:,labelidx) == 2);
yposg = ypos(datafilt(:,labelidx) == 2);
zposg = zpos(datafilt(:,labelidx) == 2);
clear xpos ypos zpos
Pr = [xposr yposr zposr];
clear xposr yposr zposr;
Pg = [xposg yposg zposg];
clear xposg yposg zposg;

% find pair wise distances
disp('Calculating pair-wise distances...');
if direction == 1
    pdistaa = pdist(Pr);
    pdistab = pdist2(Pr,Pg);
elseif direction == 2
    pdistaa = pdist(Pg);
    pdistab = pdist2(Pg,Pr);
end
pdistaa = squareform(pdistaa);
pdistaa(pdistaa == 0) = NaN;
%pdistab(pdistab == 0) = NaN;

%count # neighbors vs. r 
Naa = zeros(size(pdistaa,1),length(r));
Nab = zeros(size(pdistab,1),length(r));
disp('counting neighbors...')
for a = 1:length(r)
    Naa(:,a) = sum((pdistaa <= r(a)),2);
    Nab(:,a) = sum((pdistab <= r(a)),2);
end

%normalize by max N and max R
Daa = zeros(size(pdistaa,1),length(r));
Dab = zeros(size(pdistab,1),length(r));
disp('normalizing...');
for a = 1:size(pdistaa,1)
    Daa(a,:) = Naa(a,:)./Naa(a,end).*(rmax./r).^2;
    Dab(a,:) = Nab(a,:)./Nab(a,end).*(rmax./r).^2;
end
%Daa(isnan(Daa)) = 0;
%Dab(isnan(Dab)) = 0;

%perform spearman rank correlation on each point
disp('Finding spearman rank correlation coefficients...')
rho = zeros(size(pdistaa,1),1);

for b = 1:size(rho,1)
    disp(b)
    rho(b) = corr(Daa(b,:)',Dab(b,:)','Type','Spearman');
end
rho(isnan(rho)) = 0;

%calculate colocalization coeff for each point
minpdistab = min(pdistab,[],2);
C = rho.*exp(-1*minpdistab./rmax);

%show visualization
map = colormap; 
%close;
minC = min(C);
maxC = max(C);
lutr = interp1(linspace(minC,maxC,length(map)),map(:,1),C);
lutg = interp1(linspace(minC,maxC,length(map)),map(:,2),C);
lutb = interp1(linspace(minC,maxC,length(map)),map(:,3),C);
lut = [lutr lutg lutb];
set(handles.figure1,'CurrentAxes',handles.axes1);
cla
%handles.axes1.ColorOrder = lut; hold on
if direction == 1
    %for a = 1:10
    for a = 1:length(Pr)
        disp(a);
        plot3(Pr(a,1),Pr(a,2),Pr(a,3),'Marker','.','Color',lut(a,:)); hold on;
    end
    hold off;
elseif direction == 2
    for a = 1:length(Pg)
        disp(a);
        plot3(Pg(a,1),Pg(a,2),Pg(a,3),'Marker','.','Color',lut(a,:)); hold on;
    end
    hold off;
end

axisequal = get(handles.axisequal,'value');
if axisequal
    axis('equal');
end
c = colorbar;
currlabels = c.TickLabels;
numclabels = length(currlabels);
cinc = (max(C)-min(C))/(numclabels-1);
clabels = minC:cinc:maxC;
clabels2 = arrayfun(@num2str,clabels,'uniformoutput',0);
c.TickLabels = clabels2;

datafilt(datafilt(:,labelidx) == direction,colocidx) = C;
datamat((datamat(:,labelidx) == direction & keepindx'),colocidx) = C;
setappdata(handles.filter,'datafilt',datafilt);
setappdata(handles.filter,'datamat',datamat);
datatab = handles.datatable.Data;
datatab(colocidx,:) = [min(datafilt(:,colocidx)) max(datafilt(:,colocidx))];
handles.datatable.Data = datatab;


% --- Executes on button press in fitsurf.
function fitsurf_Callback(hObject, eventdata, handles)
% hObject    handle to fitsurf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in ripley.
function ripley_Callback(hObject, eventdata, handles)
% hObject    handle to ripley (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in clusteranalysis.
function clusteranalysis_Callback(hObject, eventdata, handles)
% hObject    handle to clusteranalysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ipalmcluster;


% --- Executes on button press in xyz.
function xyz_Callback(hObject, eventdata, handles)
% hObject    handle to xyz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
paramlist = getappdata(handles.loaddata,'paramlist');
xposidx = find(contains(paramlist,'X_Position'),1);
yposidx = find(contains(paramlist,'Y_Position'),1);
zposidx = find(contains(paramlist,'Unwrapped_Z'),1);
handles.xvals.Value = xposidx;
handles.yvals.Value = yposidx;
handles.zvals.Value = zposidx+1;


% --- Executes on button press in groupxyz.
function groupxyz_Callback(hObject, eventdata, handles)
% hObject    handle to groupxyz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
paramlist = getappdata(handles.loaddata,'paramlist');
xposidx = find(contains(paramlist,'Group_X_Position'),1);
yposidx = find(contains(paramlist,'Group_Y_Position'),1);
zposidx = find(contains(paramlist,'Unwrapped_Group_Z'),1);
handles.xvals.Value = xposidx;
handles.yvals.Value = yposidx;
handles.zvals.Value = zposidx+1;


% --- Executes on button press in nearestneighbor.
function nearestneighbor_Callback(hObject, eventdata, handles)
% hObject    handle to nearestneighbor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
datafilt = getappdata(handles.filter,'datafilt');
datamat = getappdata(handles.filter,'datamat');
keepindx = getappdata(handles.filter,'keepindx');
paramlist = getappdata(handles.loaddata,'paramlist');
xposidx = find(contains(paramlist,'X_Position'),1);
yposidx = find(contains(paramlist,'Y_Position'),1);
zposidx = find(contains(paramlist,'Unwrapped_Z'),1);
labelidx = find(contains(paramlist,'Label_Set'),1);
xpos = datafilt(:,xposidx);
ypos = datafilt(:,yposidx);
zpos = datafilt(:,zposidx);
label = datafilt(:,labelidx);
labels = unique(label);
output = inputdlg({'Distance from label','Distance to label'},'Choose which labels for NN Distance');
distfromlab = str2double(output{1});
disttolab = str2double(output{2});
if ~ismember(distfromlab,labels)
    error('Distance from label not present');
end
if ~ismember(disttolab,labels)
    error('Distance to label not present');
end

if distfromlab == disttolab
    xyz = [xpos(label == distfromlab) ypos(label == distfromlab) zpos(label == distfromlab)];
    distmat = pdist(xyz);
    distmat = squareform(distmat);
elseif distfromlab ~= disttolab
    xyz1 = [xpos(label == distfromlab) ypos(label == distfromlab) zpos(label == distfromlab)];
    xyz2 = [xpos(label == disttolab) ypos(label == disttolab) zpos(label == disttolab)];
    distmat = pdist2(xyz1,xyz2);
    %distmat = sqrt(distmat);
end
nndist = min(distmat,[],2);
NNidx = find(contains(paramlist,'NN_Distance'),1);
datafilt(label == distfromlab,NNidx) = nndist;
labelmat = datamat(:,labelidx);
datamat((labelmat == distfromlab) & keepindx',NNidx) = nndist;
setappdata(handles.filter,'datafilt',datafilt);
setappdata(handles.filter,'datamat',datamat);
datatab = handles.datatable.Data;
datatab(NNidx,:) = [min(datafilt(:,NNidx)) max(datafilt(:,NNidx))];
handles.datatable.Data = datatab;


% --- Executes on button press in rotate.
function rotate_Callback(hObject, eventdata, handles)
% hObject    handle to rotate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
datafilt = getappdata(handles.filter,'datafilt');
paramlist = getappdata(handles.loaddata,'paramlist');
output = inputdlg('Rotate Angle (degrees)','Rotate Filtered Data');
thet = str2double(output{1});
xposidx = find(contains(paramlist,'X_Position'),1);
yposidx = find(contains(paramlist,'Y_Position'),1);
xpos = datafilt(:,xposidx);
ypos = datafilt(:,yposidx);
p1 = [xpos ypos];
xcenter = (max(xpos) - min(xpos))/2;
ycenter = (max(ypos) - min(ypos))/2;
center = repmat([xcenter ycenter],length(xcenter),1);
p2 = p1-center;
R = [cosd(thet) -1*sind(thet); sind(thet) cosd(thet)];
p3 = R*p2';
p4 = p3'+center;
datafilt(:,xposidx) = p4(:,1);
datafilt(:,yposidx) = p4(:,2);
setappdata(handles.filter,'datafilt',datafilt);


% --- Executes on button press in revz.
function revz_Callback(hObject, eventdata, handles)
% hObject    handle to revz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of revz


% --- Executes on button press in removepts.
function removepts_Callback(hObject, eventdata, handles)
% hObject    handle to removepts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
xvalnum = handles.xvals.Value;
yvalnum = handles.yvals.Value;
paramlist = getappdata(handles.loaddata,'paramlist');
if ~strcmp(paramlist{xvalnum},'X_Position') && ~strcmp(paramlist{yvalnum},'Y_Position')
    error('data should be plotted in xyz space')
end
keepindx = getappdata(handles.filter,'keepindx');
datamat = getappdata(handles.filter,'datamat');
axes(handles.axes1);
revz = handles.revz.Value;
if ~revz
    view([0 90])
elseif revz
    view([0 -90])
end
rect = getrect;
xvalmin = rect(1); 
yvalmin = rect(2);
xvalmax = rect(1)+rect(3);
yvalmax = rect(2)+rect(4);
indxthrow = datamat(:,xvalnum)>xvalmin & datamat(:,xvalnum)< xvalmax & datamat(:,yvalnum)>yvalmin & datamat(:,yvalnum)<yvalmax;
keepindx(indxthrow) = 0;
datafilt = datamat(keepindx,:);
setappdata(handles.filter,'keepindx',keepindx);
setappdata(handles.filter,'datafilt',datafilt);
plotvals_Callback(handles.plotvals, eventdata, handles);
if ~revz
    view([0 90])
elseif revz
    view([0 -90])
end







