figh=tiledlayout(1,2);figh.TileSpacing = 'none';
%% range plot example
ax1=nexttile;
series={'Réglisse','Chocolat','Caramel','Nougat','Meringue','Guimauve'};
values=[ 2.4 5 3; 1.2 4.3 3; 2 5.25 3; 4 6 5; 3.2 6.4 5; 48 550 100];
options.patchHeigth=0.3;options.YStep=3;options.formatstr={'%1.1f';'%0.2f';'%1.2f';'%1.0f';'%1.1f';'%1.0f'};
h=rangedotplot(ax1,series,values,options);title('Range Plot');
%% sunburst plot example
ax2=nexttile;
[patchhl, labelhl]=sunburstplot('testdata.csv',[],"Sales",2,1,0.7,1,340);title(['Sunburst or' sprintf('\n') 'Polar Treemap Plot']);
%sunburstplot('testdata.xlsx',{'r','g','b'},"Books",1);
%sunburstplot('testdata.xlsx',{[1 0 0],[0 1 0],[1 1 0]},"Sports",[],1);
%[patchhl, labelhl]=sunburstplot('testdata.xlsx',[],"Coffee",0,[],[],0,192); %takes sometime if plotting all four columns of coffee dataset
set(gcf, 'Color', 'w');