function [dsc] = Plot2PhaseSegCompLength (folder, RunID, xphs, ind)
%
% dsc = PlotSegCompLength (folder, RunID, ind)
%
% example: [dsc] = Plot2PhaseSegCompLength (folder, RunID, [1,2]);
%
% plots the segregation-compaction length of a two-phase system with the
% range of phase fractions from simulation marked with grey rectangle
% 
% INPUTS
% folder    folder name where runs are stored
% RunID     name of run
% xphs      which phase you want on the x axis [scalar]
% ind       index of the seg-comp length that you want to plot [2x1]
%
% OUTPUTS
% dsc       segregation-compaction lengths
%
% YQW, 26 May 2021

% check inputs
if nargin<3, xphs = 2; end

% extract file names
[fn, fp] = GetOutputMatFiles(folder, RunID);

load(fp, 'eta0','d0','A','B','C');
load(fn{end}, 'f');

f1 = linspace(0,1,1001);
fdsc  = [f1; 1-f1];

dsc = SegCompLength(fdsc, eta0, d0, A, B, C);

figure;

% plot general segregation-compaction length relationship
if nargin>4
    semilogy(fdsc(xphs,:), squeeze(dsc(ind(1),ind(2),:)));
else
    semilogy(fdsc(xphs,:), squeeze(max(dsc,[],[1,2])));
end

% plot the range of phase fractions from simulation
hold on;
ylimits = ylim;
frange  = [min(f(xphs,:,:),[],'all'), max(f(xphs,:,:),[],'all')];
r = rectangle('Position',[frange(1), ylimits(1), diff(frange), diff(ylimits)],...
    'FaceColor', 0.8*ones(1,3), 'EdgeColor', 'none');
uistack(r,'bottom');
hold off;

xlabel('f2'); ylabel('seg-comp length [m]');




end