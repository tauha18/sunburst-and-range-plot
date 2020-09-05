function hl = sunburstplot(filename,varargin)
% Author: Muhammad Tauha Ali (August 2020)
%% sunburst:
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
%  varargin{5}: Root circle radius (change to accomodate Root label if needed) (default=0.5)
%  other options that can be changed in code (without detailed understanding) are:
%  1. Sector width (default=1) 2. Space between rings (default=0)
%  3. Starting angle in radians for first sector (default 0) 4. Ring transparency range (default: [0.5 0.2])
%  5. Font name (default: Arial) 6. Font size (default: 11)
%
% Arguments: (output)
%  hl - handle to the sector patches. Each ring patches have their own
%  dimension i.e. First dimension represent first ring, second dimension
%  represent second ring and so on
% Examples:
%   sunburstplot(filename,{'w','r',g','b','m'},"Books");
%   patchhl=sunburstplot(filename,{[ .945 .345 .329],[ .376 .741 .408]},'Books',1);

% Default values if no variable arguments are provided
clrmp=[];
Rootlabel="";%Label to appear at root
valflag=2;% value display flag
rotflag=1;% label rotation flag
rootrad=0.5;%Radius of root label circle
secwidth=1.25;%width of circle rings
ringgap=0.0;%space between subsequent rings
stangle=pi/2*0;%starting angle of first sector
alpharange=[0.5 0.2]; %patch transparency range
fontnm='Arial'; %text font name
fontpt=11; %font size

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
for i=1:size(datatree,2) %convert cells to string arrays
    if ~isnumeric(eval(['datatree.Var' num2str(i)]))
        eval(['datatree.Var' num2str(i) '=string(datatree.Var' num2str(i) ');']);
    else
        if i==size(datatree,2)
            valcolflag=1; %last column contain numbers
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
        case 3
            if ~isempty(varargin{3})
                valflag = varargin{3};
            end
        case 4
            if ~isempty(varargin{4})
                rotflag = varargin{4};
            end
        case 5
            if ~isempty(varargin{5})
                rootrad = varargin{5};
            end
        otherwise
            error('Too many input arguments');
    end
end
%if last column does not contain numbers then add column of 1 for counting labels
if ~exist('valcolflag','var')
    datatree(:,size(datatree,2)+1)=num2cell(ones(size(datatree,1),1));
end
colno=size(datatree,2)-1; %number of rings

% collect sector labels and remove any empty labels
seclabels=cell(1,colno);
for i=1:colno
    seclabels(i)={unique(string(table2cell(datatree(:,i))))};
    seclabels{i}(seclabels{:,i}=="" | seclabels{:,i}==" " | ismissing(seclabels{:,i}))=[];
end

%specify sector order if desired otherwise ascending label-wise
% seclabels{1}=seclabels{1}([1 2 4 3],1);
% seclabels{2}=seclabels{2}([2 1 3],1);
% seclabels{3}=seclabels{3}([2 1],1);

branches = size(seclabels{1},1); % number of branches in first ring
% generate or correct branches colormap
if isempty(clrmp)
    clrmp=num2cell(distinguishable_colors(branches),2);
elseif ~(size(clrmp,2)==branches) %ignore colormap if supplied colormap does not have each branch color
    clrmp=num2cell(distinguishable_colors(branches),2); %color is determined  by number of branches
    warning('Supplied colormap is ignored because of incorrect color values. Supply color value for each branch only');
end
hl=zeros(branches,1);%patch handles variable initialization
patchalphab=max(alpharange); %transperancy will increase outward
for i=1:colno
    switch i
        case 2
            leaves1 = size(seclabels{2},1);
            hl=zeros(branches,leaves1+1);%patch handles variable initialization
            patchalphal1=min(alpharange);
        case 3
            leaves2 = size(seclabels{3},1);
            hl=zeros(branches,leaves1+1,leaves2+1);%patch handles variable initialization
            patchalphal1=mean(alpharange);patchalphal2=min(alpharange);
        case 4
            leaves3 = size(seclabels{4},1);
            hl=zeros(branches,leaves1+1,leaves2+1,leaves3+1);%patch handles variable initialization
            dum=linspace(max(alpharange),min(alpharange),4);
            patchalphal1=dum(2);patchalphal2=dum(3);patchalphal3=dum(4);
    end
end

% sunburst plotting code
fractang=stangle;%starting angle for first sector in each ring
for i = 1:branches %sector plotting is done branch-wise
    subdata=datatree(string(table2cell(datatree(:,1)))==seclabels{1}(i),:);
    frac=nansum(table2array(subdata(:,end)))./nansum(table2array(datatree(:,end)));
    fractang = [fractang,fractang+frac.*(2*pi)];
    
    r0 = rootrad;r1 = rootrad+(secwidth-ringgap);%width of sector
    cl = clrmp{i};% sector color-color remains same for branches but transparency change
    
    if (fractang(2)-fractang(1))~=0 %plot if sector exists
        hl(i,1)=polsect(fractang(1),fractang(2),r0,r1,cl,patchalphab); %plot sector
        labelRadius=(r0+r1)/2;centerTheta=mean(fractang);
        [halign,valign,rotang]=getAlignmentFromAngle(centerTheta,rotflag,45);%get text rotation and text alignment
        [xtext,ytext] = pol2cart(centerTheta,labelRadius);
        switch valflag
            case 0
                labelf=seclabels{1}(i);
            case 1
                labelf=[seclabels{1}(i) num2str(round(nansum(table2array(subdata(:,end)))))];
            case 2
                labelf=[seclabels{1}(i) num2str(round(frac*100))+"%"];
        end
        text(xtext,ytext,labelf,'fontname',fontnm,'fontsize',fontpt,'rotation',rotang,'HorizontalAlignment',halign,'VerticalAlignment',valign);
    end
    if exist('leaves1','var') %if level1 leaves column is present
        fractangl1=fractang(1);phl1id=2;
        for j=1:leaves1
            r0 = rootrad+secwidth;r1 = rootrad+secwidth+(secwidth-ringgap);%width of sector
            subdata1=subdata(string(table2cell(subdata(:,2)))==seclabels{2}(j),:);
            frac=nansum(table2array(subdata1(:,end)))./nansum(table2array(subdata(:,end)));frac(isnan(frac))=0;
            fractangl1 = [fractangl1,frac.*(fractang(2)-fractang(1))+fractangl1];
            
            if (fractangl1(2)-fractangl1(1))~=0 %plot if sector exists
                hl(i,phl1id,1)=polsect(fractangl1(1),fractangl1(2),r0,r1,cl,patchalphal1); %plot sector
                labelRadius=(r0+r1)/2;centerTheta=mean(fractangl1);
                [halign,valign,rotang]=getAlignmentFromAngle(centerTheta,rotflag,22.5);
                [xtext,ytext] = pol2cart(centerTheta,labelRadius);
                switch valflag
                    case 0
                        labelf=seclabels{2}(j);
                    case 1
                        labelf=[seclabels{2}(j) num2str(round(nansum(table2array(subdata1(:,end)))))];
                    case 2
                        labelf=[seclabels{2}(j) num2str(round(frac*100))+"%"];
                end
                text(xtext,ytext,labelf,'fontname',fontnm,'fontsize',fontpt,'rotation',rotang,'HorizontalAlignment',halign,'VerticalAlignment',valign);
            else
                hl(i,phl1id,1)=0;
            end
            if exist('leaves2','var') %if Level 2 leaves column exist
                fractangl2=fractangl1(1);phl2id=2;
                for k=1:leaves2
                    r0 = rootrad+secwidth*2;r1 = rootrad+secwidth*2+(secwidth-ringgap);%width of sector
                    subdata2=subdata1(string(table2cell(subdata1(:,3)))==seclabels{3}(k),:);
                    frac=nansum(table2array(subdata2(:,end)))./nansum(table2array(subdata1(:,end)));frac(isnan(frac))=0;
                    fractangl2 = [fractangl2,frac.*(fractangl1(2)-fractangl1(1))+fractangl2];
                    
                    if (fractangl2(2)-fractangl2(1))~=0 %plot if sector exists
                        hl(i,phl1id,phl2id)=polsect(fractangl2(1),fractangl2(2),r0,r1,cl,patchalphal2); %plot sector
                        labelRadius=(r0+r1)/2;centerTheta=mean(fractangl2);[halign,valign,rotang]=getAlignmentFromAngle(centerTheta,rotflag,22.5/2);
                        [xtext,ytext] = pol2cart(centerTheta,labelRadius);
                        switch valflag
                            case 0
                                labelf=seclabels{3}(k);
                            case 1
                                labelf=[seclabels{3}(k) num2str(round(nansum(table2array(subdata2(:,end)))))];
                            case 2
                                labelf=[seclabels{3}(k) num2str(round(frac*100))+"%"];
                        end
                        text(xtext,ytext,labelf,'fontname',fontnm,'fontsize',fontpt,'rotation',rotang,'HorizontalAlignment',halign,'VerticalAlignment',valign);
                    else
                        hl(i,phl1id,phl2id)=0;
                    end
                    if exist('leaves3','var') %if Level 3 leaves column exist
                        fractangl3=fractangl2(1);phl3id=2;
                        for m=1:leaves3
                            r0 = rootrad+secwidth*3;r1 = rootrad+secwidth*3+(secwidth-ringgap);%width of sector
                            subdata3=subdata2(string(table2cell(subdata2(:,4)))==seclabels{4}(m),:);
                            frac=nansum(table2array(subdata3(:,end)))./nansum(table2array(subdata2(:,end)));frac(isnan(frac))=0;
                            fractangl3 = [fractangl3,frac.*(fractangl2(2)-fractangl2(1))+fractangl3];
                            
                            if (fractangl3(2)-fractangl3(1))~=0 %plot if sector exists
                                hl(i,phl1id,phl2id,phl3id)=polsect(fractangl3(1),fractangl3(2),r0,r1,cl,patchalphal3); %plot sector
                                labelRadius=(r0+r1)/2;centerTheta=mean(fractangl3);
                                %[halign,valign,rotang]=getAlignmentFromAngle(centerTheta,rotflag,22.5/2/2);
                                halign='center';valign='middle';angle=mod(round(centerTheta*180/pi/(22.5/4))*22.5/4,360);
                                if angle>90 && angle<270
                                    rotang=centerTheta*180/pi-180;
                                else
                                    rotang=centerTheta*180/pi;
                                end
                                [xtext,ytext] = pol2cart(centerTheta,labelRadius);
                                switch valflag
                                    case 0
                                        labelf=seclabels{4}(m);
                                    case 1
                                        labelf=[seclabels{4}(m) num2str(round(nansum(table2array(subdata3(:,end)))))];
                                    case 2
                                        labelf=[seclabels{4}(m) num2str(round(frac*100))+"%"];
                                end
                                text(xtext,ytext,labelf,'fontname',fontnm,'fontsize',fontpt,'rotation',rotang,'HorizontalAlignment',halign,'VerticalAlignment',valign);
                            else
                                hl(i,phl1id,phl2id,phl3id)=0;
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
text(0,0,Rootlabel,'fontname',fontnm,'fontsize',fontpt+2,'fontweight','Bold','HorizontalAlignment','center','VerticalAlignment','middle');
axis equal;ax=gca;ax.XAxis.Visible='off';ax.YAxis.Visible='off';
%{
% if i==rings
%     legend1 = legend(legtext);
%     wi = legend1.Position(3);
%     Xlm = xlim;
%     widx = diff(Xlm);
%     unitwi = widx.*wi;
%     xlim([Xlm(1),Xlm(2)+unitwi])
% end
%}
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
angle = (180/pi)*angle;

% Round the angles to the nearest modangle in degrees
angle = mod(round(angle/modangle)*modangle,360);

% Determine the horizontal alignment (center works best usually)
if angle == 90 || angle == 270
    halign = 'center';
elseif angle > 90 && angle < 270
    halign = 'center';%'right';
else
    halign = 'center';%'left';
end

% Determine the vertical alignment (middle works best usually)
if angle == 0 || angle == 180
    valign = 'middle';
elseif angle > 180 && angle < 360
    valign = 'middle';%'top'
else
    valign = 'middle';%'bottom'
end

% Calculate rotation angle
if angle>=0 && angle<180
    rotangle=angle-90;
elseif angle>=180 && angle<270
    rotangle=angle-180;
elseif angle>=270 && angle<360
    rotangle=angle-270;
end
if ~rotflag
    rotangle=0;
end
end