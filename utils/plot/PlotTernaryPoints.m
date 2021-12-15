function [] = PlotTernaryPoints (f, v, MkSize)
%
% PlotTernaryPoints (f, v, MkSize, fName, vName)
% 
% plot points on a ternary diagram for a three phase system (NPHS = 3)
% 
% INPUTS
% f         phase fraction of points [NPHS x Npts]
% v         variable value to plot [Npts x 1]
% MkSize    marker size [Npts x 1 or 1]

if nargin<3, MkSize = 100; end

% check size of f
if size(f,1)==3 && size(f,2)~=3
    f = f';
end

[X,Y]  =  terncoords(f(1,:),f(2,:),f(3,:));

hold on;

a = scatter(X, Y, MkSize, v, 'filled');
hold on
uistack(a, 'top');


end