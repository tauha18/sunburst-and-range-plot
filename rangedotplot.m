
function varargout = rangedotplot(axh,series,values,varargin)

%% rangedotplot:
% create a figure with patch bars to represent data series with value range
% and optionally specify a value within the range to show as dot (added by tauha)
%
% usage: rangedotplot(axeshandle,series,values)
% usage: handleobj = rangedotplot(axeshandle,series,values)
% usage: handleobj = rangedotplot(axeshandle,series,values,options)
%
% Arguments: (input)
%  series - cell array - contains data series names.
%  values - numeric array - contains min and max values for each data serie.
%  options - (OPTIONAL) structure - contains options.  Valid properties:
%       'patchHeigth', 'displayValues', 'YStep',
%       'linecolor','formatstr' (added options by tauha)
%
%  'patchHeigth' - heigth of the patch (between 0 and 1...) (default=0.2)
%  'displayValues' - true --> display the values (default)
%                     false --> or not
%  'YStep' - heigth of the step between two patches (default=1)
%  'linecolor' - color of each range patch
%  (default uses distinguishable_colors code to generate colors)
%  'formatstr' - format of numeric values for each range (default=%4.2f)
%
% Arguments: (output)
%  handleFig - handle to the figure and handle to patch objects
%
%  If no output argument is provided for, then the figure
%  will only be displayed.
%
% Example usages:
%  series={'Réglisse','Chocolat','Caramel','Nougat','Meringue','Guimauve'};
%  values=[ 2.4 5 ; 1.2 4.3 ; 2 5.25 ; 4 6; 3.2 6.4; 0.8 5.5];
%  options.patchHeigth=0.6;
%  rangedotplot_mod(series,values,options)
%
% Copyright Fabien Baillon 29/07/2017 for plotStickers -modified by Tauha 17/08/2020
% Fabien.Baillon@mines-albi.fr

nSeries=numel(series);
% defaults for all optional parameters
options.linecolor =distinguishable_colors(nSeries);
options.formatstr=repmat({'%4.2f'},nSeries,1);
options.patchHeigth=0.2;
options.displayValues=true;
options.YStep=1;

% User can override default options
if(nargin>2)
    userOptions=varargin{1};
    if(isstruct(userOptions))
        varF=fieldnames(userOptions);
        for iF=1:numel(varF)
            if(isfield(options,varF{iF}))
                options.(varF{iF})=userOptions.(varF{iF});
            end
        end
    end
end
normvalue=(values(:,1)+values(:,2))/2;values=values./normvalue;%--modified code--
% Data preparation for the patch function
pX=[values(:,1:2) flip(values(:,1:2),2)];

stepY=options.YStep+5*(options.patchHeigth-options.YStep).*(options.patchHeigth>options.YStep);
stepsY=stepY*[1:nSeries]';
pY=repmat(stepsY,1,4)+repmat(options.patchHeigth*[-1 -1 1 1]/2,nSeries,1);

% Creation of the figure
box(axh,'off');
axYlabels=cell(nSeries+2,1);
axYlabels(2:nSeries+1)=series(:);
set(axh,'Xlim',[min(min(values))-0.2,max(max(values))+0.3],'YLim',[stepY-options.patchHeigth,stepY*nSeries+options.YStep],'YTick',[],'YTickLabel',[]);
axh.YAxis.Visible = 'off';axh.XAxis.Visible = 'off';

% Adding data as patch stickers
hold on
for iL=1:nSeries
    h(iL)=patch(pX(iL,:),pY(iL,:),options.linecolor(iL,:),'linestyle','none');
    %--modified code--
    if size(values,2)==3 %if values contain the optional third column for dot plotting
        plot(values(iL,3),mean(pY(iL,[1 3])),'o','color',options.linecolor(iL,:),'MarkerFaceColor',options.linecolor(iL,:),'MarkerSize',8);
    end
end
%--modified code--
set(h(:),'Facealpha',0.25);
% Diplaying the values
if(options.displayValues)
    for i=1:nSeries
        leftLabels(i,1)={sprintf([options.formatstr{i} '   '],values(i,1).*normvalue(i))};
        rightLabels(i,1)={sprintf(['   ' options.formatstr{i}],values(i,2).*normvalue(i))};
        if size(values,2)==3
            dotLabels(i,1)={sprintf(['  ' options.formatstr{i}],values(i,3).*normvalue(i))};
        end
    end
    text(values(:,2), stepsY(:), rightLabels(:));
    text(values(:,1), stepsY(:), leftLabels(:), 'HorizontalAlignment','right');
    if size(values,2)==3
        text(values(:,3), stepsY(:), dotLabels(:));
    end
end

% Adding labels for each bar
text((values(:,1)+values(:,2))/2,stepsY(:)+1,series,'HorizontalAlignment','center');
%set(axh,'YTick',stepY*(0:nSeries+1),'YTickLabel', axYlabels);
hold off
%-------------------
if(nargout)
    varargout{1}=h';
end

end