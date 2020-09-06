function [phl, thl] = sunburstplot(filename,varargin)
% Author: Muhammad Tauha Ali (August 2020)
%% Sunburst Plot/Polar Treemap Plot:
% datatree: cell type data. Each column denotes ring level. (do not include header)
% Column format: 1. branches 2. Leaves (only three levels allowed at this time)
% Every column should contain sector labels as strings
% if last column contain numbers then they are used to compute sector angles
%
% Arguments: (input)
%  varargin{1}: Cell of colors as alphabets or one row of 3 RGB values for each branch in datatree
%  (default uses distinguishable_colors code to generate colors)
%  varargin{2}: Root label (default: "")
%  varargin{3}: Value display flag: 0. No values 1. numbers 2. percentages (default)
%  varargin{4}: Label rotation flag: 0. do not rotate 1. rotate label according to sector angles (default)
%  2. rotate labels so that sector labels appear radially (Labels in 2nd leaves level and above are radially plotted by default)
%  varargin{5}: Root circle radius (change to accomodate Root label if needed) (default=0.5)
%  other options present in the start that can be changed (without detailed code understanding) are:
%  varargin{6}: Plotting direction: 0. clockwise 1. couter clockwise (default = 1)
%  varargin{7}: Starting angle in radians for first sector (default 0)
%  1. Sector width (default=1) 2. Space between rings (default=0)
%  3. Ring color transparency range (default: [0.5 0.2])
%  4. Font name (default: Arial) 5. Font size (default: 11)
%
% Arguments: (output)
%  phl - handle to the sector patches. Each ring patches have their own
%  dimension i.e. First dimension represent first ring, second dimension
%  represent second ring and so on
%  thl - handle to the sector labels in similar manner to patch handles
% Examples:
%   sunburstplot(filename,{'w','r',g','b','m'},"Books");
%   [patchhl labelhl]=sunburstplot(filename,{[ .945 .345 .329],[ .376 .741 .408]},'Books',1);

% Default values if no variable arguments are provided
clrmp=[];
Rootlabel="";   % Label to appear at root
valflag=2;      % value display flag
rotflag=1;      % label rotation flag
rootrad=0.5;    % Radius of root label circle
stdir=1;        % patch plotting direction
stangle=0*pi/180; % starting angle of first sector

secwidth=1.25;  % width of circle rings
ringgap=0.0;    % space between subsequent rings
alpharange=[0.5 0.2]; % patch transparency range
fontnm='Arial'; %text font name
fontpt=11;      % font size

% data reading
if strcmp(filename(end-3:end),'xlsx') || strcmp(filename(end-2:end),'xls')
    [~,~,datatree]=xlsread(filename,-1);datatree(cellfun(@(x) isequaln(x,NaN),datatree))={''};
    Varlabel=cell(1,size(datatree,2));
    for i=1:size(datatree,2)
        Varlabel(i)={['Var' num2str(i)]};
    end
    datatree=cell2table(datatree,'VariableNames',Varlabel);
else
    datatree=readtable(filename);
end
for i=1:size(datatree,2) % convert cells to string arrays
    if ~isnumeric(eval(['datatree.Var' num2str(i)]))
        eval(['datatree.Var' num2str(i) '=string(datatree.Var' num2str(i) ');']);
    else
        if i==size(datatree,2)
            valcolflag=1; % last column contain numbers
        else
            error('All Columns should contain labels except last column which can contain numbers');
        end
    end
end

% varargin assignment
for i=1:length(varargin)
    switch i
        case 1
            if ~isempty(varargin{1})
                clrmp = varargin{1};
            end            
        case 2
            if ~isempty(varargin{2})
                Rootlabel = varargin{2};
            end
            if ~(isa(Rootlabel,'string') || isa(Rootlabel,'char'))
                error('Label should either string or character type');
            end
        case 3
            if ~isempty(varargin{3})                
                valflag = varargin{3};
                if ~ismember(valflag,[0 1 2])
                    error('Label value flag should be either 0, 1 or 2');
                end
            end
        case 4
            if ~isempty(varargin{4})
                rotflag = varargin{4};
            end
            if ~ismember(rotflag,[0 1 2])
                error('Label rotation flag should be either 0, 1 or 2');
            end
        case 5
            if ~isempty(varargin{5})
                rootrad = varargin{5};
            end
            if ~isnumeric(rootrad) || rootrad<0
                error('Root radius should be a positive number');
            end
        case 6            
            if ~isempty(varargin{6})
                stdir = varargin{6};
            end
            if ~ismember(stdir,[0 1])
                error('Plotting direction flag should be either 0 or 1');
            end
        case 7
            if ~isempty(varargin{7})
                if varargin{7}<0 || varargin{7}>360
                    error('Starting angle should be between 0 and 360 degrees');
                else
                stangle = varargin{7}*pi/180;
                end
            end
        otherwise
            error('Too many input arguments');
    end
end
% if last column does not contain numbers then add column of 1 for counting labels
if ~exist('valcolflag','var')
    datatree(:,size(datatree,2)+1)=num2cell(ones(size(datatree,1),1));
end
colno=size(datatree,2)-1; % number of rings

% collect sector labels and remove any empty labels
seclabels=cell(1,colno);
for i=1:colno
    seclabels(i)={unique(string(table2cell(datatree(:,i))))};
    seclabels{i}(seclabels{:,i}=="" | seclabels{:,i}==" " | ismissing(seclabels{:,i}))=[];
    if ~stdir
        seclabels{i}=seclabels{i}(fliplr(1:length(seclabels{i})),1);
    end
end

% specify sector order if desired otherwise ascending label-wise
% seclabels{1}=seclabels{1}([1 2 4 3],1);
% seclabels{2}=seclabels{2}([2 1 3],1);
% seclabels{3}=seclabels{3}([2 1],1);

branches = size(seclabels{1},1); % number of branches in first ring
% generate or correct branches colormap
if isempty(clrmp)
    clrmp=num2cell(distinguishable_colors(branches),2);
elseif ~(size(clrmp,2)==branches) % ignore colormap if supplied colormap does not have each branch color
    clrmp=num2cell(distinguishable_colors(branches),2); % color is determined  by number of branches
    warning('Supplied colormap is ignored because of incorrect color values. Supply color value for each branch only');
end
phl=zeros(branches,1);thl=zeros(branches,1);% patch and label handles variable initialization
patchalphab=max(alpharange); % transperancy will increase outward
for i=1:colno
    switch i
        case 2
            leaves1 = size(seclabels{2},1);
            phl=zeros(branches,leaves1+1);thl=zeros(branches,leaves1+1);% patch and label handles variable initialization
            patchalphal1=min(alpharange); % patch transperancy assignment
        case 3
            leaves2 = size(seclabels{3},1);
            % patch and label handles variable initialization
            phl=zeros(branches,leaves1+1,leaves2+1);thl=zeros(branches,leaves1+1,leaves2+1);
            % patch transperancy assignment
            patchalphal1=mean(alpharange);patchalphal2=min(alpharange);
        case 4
            leaves3 = size(seclabels{4},1);
            % patch and label handles variable initialization
            phl=zeros(branches,leaves1+1,leaves2+1,leaves3+1);thl=zeros(branches,leaves1+1,leaves2+1,leaves3+1);
            % patch transperancy assignment
            dum=linspace(max(alpharange),min(alpharange),4);
            patchalphal1=dum(2);patchalphal2=dum(3);patchalphal3=dum(4);
    end
end

% sunburst plotting code
fractang=stangle;% starting angle for first sector in each ring
for i = 1:branches % sector plotting is done branch-wise
    subdata=datatree(string(table2cell(datatree(:,1)))==seclabels{1}(i),:);
    frac=nansum(table2array(subdata(:,end)))./nansum(table2array(datatree(:,end)));
    fractang = [fractang,fractang+frac.*(2*pi)];
    
    r0 = rootrad;r1 = rootrad+(secwidth-ringgap);% width of sector
    cl = clrmp{i};% sector color-color remains same for branches but transparency change
    
    if (fractang(2)-fractang(1))~=0 % plot if sector exists
        phl(i,1)=polsect(fractang(1),fractang(2),r0,r1,cl,patchalphab); % plot sector
        labelRadius=(r0+r1)/2;centerTheta=mean(fractang);
        [halign,valign,rotang]=getAlignmentFromAngle(centerTheta,rotflag,45);% get text rotation and text alignment
        [xtext,ytext] = pol2cart(centerTheta,labelRadius);
        switch valflag
            case 0
                labelf=seclabels{1}(i);
            case 1
                labelf=[seclabels{1}(i) num2str(round(nansum(table2array(subdata(:,end)))))];
            case 2
                labelf=[seclabels{1}(i) num2str(round(frac*100))+"%"];
        end
        thl(i,1)=text(xtext,ytext,labelf,'fontname',fontnm,'fontsize',fontpt,'rotation',rotang,'HorizontalAlignment',halign,'VerticalAlignment',valign);
    end
    if exist('leaves1','var') % if level1 leaves column is present
        fractangl1=fractang(1);phl1id=2;
        for j=1:leaves1
            r0 = rootrad+secwidth;r1 = rootrad+secwidth+(secwidth-ringgap);% width of sector
            subdata1=subdata(string(table2cell(subdata(:,2)))==seclabels{2}(j),:);
            frac=nansum(table2array(subdata1(:,end)))./nansum(table2array(subdata(:,end)));frac(isnan(frac))=0;
            fractangl1 = [fractangl1,frac.*(fractang(2)-fractang(1))+fractangl1];
            
            if (fractangl1(2)-fractangl1(1))~=0 % plot if sector exists
                phl(i,phl1id,1)=polsect(fractangl1(1),fractangl1(2),r0,r1,cl,patchalphal1); % plot sector
                labelRadius=(r0+r1)/2;centerTheta=mean(fractangl1);
                [halign,valign,rotang]=getAlignmentFromAngle(centerTheta,rotflag,45/2);
                [xtext,ytext] = pol2cart(centerTheta,labelRadius);
                switch valflag
                    case 0
                        labelf=seclabels{2}(j);
                    case 1
                        labelf=[seclabels{2}(j) num2str(round(nansum(table2array(subdata1(:,end)))))];
                    case 2
                        labelf=[seclabels{2}(j) num2str(round(frac*100))+"%"];
                end
                thl(i,phl1id,1)=text(xtext,ytext,labelf,'fontname',fontnm,'fontsize',fontpt,'rotation',rotang,...
                    'HorizontalAlignment',halign,'VerticalAlignment',valign);
            else
                phl(i,phl1id,1)=0;thl(i,phl1id,1)=0;
            end
            if exist('leaves2','var') % if level 2 leaves column exist
                fractangl2=fractangl1(1);phl2id=2;
                for k=1:leaves2
                    r0 = rootrad+secwidth*2;r1 = rootrad+secwidth*2+(secwidth-ringgap);% width of sector
                    subdata2=subdata1(string(table2cell(subdata1(:,3)))==seclabels{3}(k),:);
                    frac=nansum(table2array(subdata2(:,end)))./nansum(table2array(subdata1(:,end)));frac(isnan(frac))=0;
                    fractangl2 = [fractangl2,frac.*(fractangl1(2)-fractangl1(1))+fractangl2];
                    
                    if (fractangl2(2)-fractangl2(1))~=0 % plot if sector exists
                        phl(i,phl1id,phl2id)=polsect(fractangl2(1),fractangl2(2),r0,r1,cl,patchalphal2); % plot sector
                        labelRadius=(r0+r1)/2;centerTheta=mean(fractangl2);
                        [halign,valign,rotang]=getAlignmentFromAngle(centerTheta,2,45/4);
                        [xtext,ytext] = pol2cart(centerTheta,labelRadius);
                        switch valflag
                            case 0
                                labelf=seclabels{3}(k);
                            case 1
                                labelf=[seclabels{3}(k) num2str(round(nansum(table2array(subdata2(:,end)))))];
                            case 2
                                labelf=[seclabels{3}(k) num2str(round(frac*100))+"%"];
                        end
                        thl(i,phl1id,phl2id)=text(xtext,ytext,labelf,'fontname',fontnm,'fontsize',fontpt,'rotation',rotang,...
                            'HorizontalAlignment',halign,'VerticalAlignment',valign);
                    else
                        phl(i,phl1id,phl2id)=0;thl(i,phl1id,phl2id)=0;
                    end
                    if exist('leaves3','var') % if level 3 leaves column exist
                        fractangl3=fractangl2(1);phl3id=2;
                        for m=1:leaves3
                            r0 = rootrad+secwidth*3;r1 = rootrad+secwidth*3+(secwidth-ringgap);% width of sector
                            subdata3=subdata2(string(table2cell(subdata2(:,4)))==seclabels{4}(m),:);
                            frac=nansum(table2array(subdata3(:,end)))./nansum(table2array(subdata2(:,end)));frac(isnan(frac))=0;
                            fractangl3 = [fractangl3,frac.*(fractangl2(2)-fractangl2(1))+fractangl3];
                            
                            if (fractangl3(2)-fractangl3(1))~=0 % plot if sector exists
                                phl(i,phl1id,phl2id,phl3id)=polsect(fractangl3(1),fractangl3(2),r0,r1,cl,patchalphal3); % plot sector
                                labelRadius=(r0+r1)/2;centerTheta=mean(fractangl3);
                                [halign,valign,rotang]=getAlignmentFromAngle(centerTheta,2,45/8);                                
                                [xtext,ytext] = pol2cart(centerTheta,labelRadius);
                                switch valflag
                                    case 0
                                        labelf=seclabels{4}(m);
                                    case 1
                                        labelf=[seclabels{4}(m) num2str(round(nansum(table2array(subdata3(:,end)))))];
                                    case 2
                                        labelf=[seclabels{4}(m) num2str(round(frac*100))+"%"];
                                end
                                thl(i,phl1id,phl2id,phl3id)=text(xtext,ytext,labelf,'fontname',fontnm,'fontsize',fontpt,...
                                    'rotation',rotang,'HorizontalAlignment',halign,'VerticalAlignment',valign);
                            else
                                phl(i,phl1id,phl2id,phl3id)=0;thl(i,phl1id,phl2id,phl3id)=0;
                            end
                            fractangl3=fractangl3(2);phl3id=phl3id+1;
                        end       
                    end
                    fractangl2=fractangl2(2);phl2id=phl2id+1;
                end
            end
            fractangl1=fractangl1(2);phl1id=phl1id+1;
        end
    end
    fractang=fractang(2);
end
text(0,0,Rootlabel,'fontname',fontnm,'fontsize',fontpt+2,'fontweight','Bold',...
    'HorizontalAlignment','center','VerticalAlignment','middle');
axis equal;ax=gca;ax.XAxis.Visible='off';ax.YAxis.Visible='off';
end
function pspatch = polsect(th0,th1,rh0,rh1,cl,patchalpha)
% This function creates a patch from polar coordinates
a1 = linspace(th0,th0);
r1 = linspace(rh0,rh1);
a2 = linspace(th0,th1);
r2 = linspace(rh1,rh1);
a3 = linspace(th1,th1);
r3 = linspace(rh1,rh0);
a4 = linspace(th1,th0);
r4 = linspace(rh0,rh0);
[X,Y]=pol2cart([a1,a2,a3,a4],[r1,r2,r3,r4]);
pspatch=patch(X,Y,cl); % Note: patch function takes text or matrix color def
set(pspatch,'FaceAlpha',patchalpha);
end

function [halign, valign,rotangle] = getAlignmentFromAngle(angle,rotflag,modangle)
% Determine the text label alignment and rotation based on the angle around the circle.

% Convert the angle to degrees
anglerad = (180/pi)*angle;

% Round the angles to the nearest modangle in degrees
anglemod = mod(round(anglerad/modangle)*modangle,360);

% Determine the horizontal alignment (center works best usually)
%{
if angle == 90 || angle == 270
    halign = 'center';
elseif angle > 90 && angle < 270
    halign = 'right';
else
    halign = 'left';
end
%}
halign='center';
% Determine the vertical alignment (middle works best usually)
%{
if angle == 0 || angle == 180
    valign = 'middle';
elseif angle > 180 && angle < 360
    valign = 'top';
else
    valign = 'bottom';
end
%}
valign='middle';

% Calculate rotation angle
switch rotflag
    case 0
        rotangle=0;
    case 1        
        if anglemod>=0 && anglemod<180
            rotangle=anglemod-90;
        elseif anglemod>=180 && anglemod<270
            rotangle=anglemod-180;
        elseif anglemod>=270 && anglemod<360
            rotangle=anglemod-270;
        end
    case 2
        if anglemod>90 && anglemod<270
            rotangle=anglerad-180;
        else
            rotangle=anglerad;
        end
end
end