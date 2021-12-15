% default plotting formats

currpath = pwd;

% load custom colormap
load('../../pantarhei/src/ocean.mat', 'ocean');

% color = [   0.0000  0.4470  0.7410 ; ...
%             0.8500  0.3250  0.0980 ; ...
%             0.9290  0.6940  0.1250 ; ...
%             0.4940  0.1840  0.5560 ; ...
%             0.0000  0.0000  0.0000];

% prepare formating options
HA = {'HorizontalAlignment','left','center','right'};
VA = {'VerticalAlignment','bottom','middle','top'};
UN = {'Units','Normalized','Centimeters'};
TX = {'Interpreter','Latex'};
TL = {'TickLabelInterpreter','Latex'};
LW = {'LineWidth',1,2,3};
FS = {'FontSize',13,16,19,24,28};
MS = {'MarkerSize',6,8,12};
LS = {'LineStyle','-','--','-.',':'};
% LC = {'Color',color};

