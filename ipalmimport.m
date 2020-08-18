function [data,paramlist,filename,pathname] = ipalmimport(pixelsize)
if nargin < 1
    pixelsize = 133;
end

[filename,pathname] = uigetfile('*.txt','Choose iPALM ASCII file');
cd(pathname);
disp('Importing iPALM ASCII data...');
data = importdata(fullfile(pathname,filename),'\t');
paramlist = data.colheaders;
for a = 1:length(paramlist)
    paramlist{a} = regexprep(paramlist{a},' ','_');
end
data = data.data;
xposcol = find(contains(paramlist,'X_Position'),1);
yposcol = find(contains(paramlist,'Y_Position'),1);
gxposcol = find(contains(paramlist,'Group_X_Position'),1);
gyposcol = find(contains(paramlist,'Group_Y_Position'),1);    
data(:,xposcol) = data(:,xposcol)*pixelsize;
data(:,yposcol) = data(:,yposcol)*pixelsize;
data(:,gxposcol) = data(:,gxposcol)*pixelsize;
data(:,gyposcol) = data(:,gyposcol)*pixelsize;    
disp('Done!');