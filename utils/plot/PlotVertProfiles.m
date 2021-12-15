function [t, xp, v, fig] = PlotVertProfiles (folder, RunID, varname, varmat, varargin)
%
% [t, xp, v, fig] = PlotVertProfiles (folder, RunID, var, varargin)
% plot vertical profiles of the model at the middle of domain or at a
% specified x point
% 
% example: PlotVertProfiles('../out/', 'olv10_plg10_bas80', {'f'})
%
% INPUTS
% folder    folder name where runs are stored
% RunID     name of run
% varname   cell vector of variable names
%               NB if you want certain axis limit options, varname has to be precise
% varmat    cell vector of variable matrices to plot
%               if mising or empty, load the variable of varname from file
%               varname, varmat must be same size
% varargin  options
%
% OUTPUTS
% t         simulation times [Nf x 1]
% xp        x positions {1, Nvar}
% v         variable profiles {1, Nvar}
%
% DEFAULT OPTIONS
% opt.fname  = '';        % extra filename info
% opt.save   = 0;         % whether to save
% 
% opt.Nplt   = 1;         % number of lines to plot
% opt.iPlt   = [];        % which timesteps to plot
% 
% opt.val    = 'mid';     % mid, mean, ind, what to plot
% opt.xdsc   = 0;         % whether to divide x by max dsc
% opt.xind   = 0;         % if val='ind', which x index to plot
% 
% opt.dfplt  = 0;         % if var == 'f', whether to plot change over time
% opt.legd   = 1;         % legend showing time?
%
% YQW, 23 Feb 2021



%  check inputs
if nargin<4, varmat = []; end


% extract file names
[fp, fn] = GetOutputMatFiles(folder, RunID);
Nf       = length(fn);

% get options
opt = defopts(Nf, varargin{:});
Nvar = length(varname);

% get profiles
[t, xp, v] = GetVertProfiles(folder, RunID, varname, varmat, opt.xind, opt.xdsc);

load(fp,'NPHS','PHS','f0');
if ~exist('PHS', 'var'), PHS = strcat({'f'}, num2str((1:NPHS)', '%d')); end

zunits = 'm';
if opt.xdsc, zunits = '\delta_{sc}'; end

% plot outputs
fig = figure;
set(gcf,'defaultaxescolororder', copper(opt.Nplt));
set(gcf,'Position', [500,100,200*NPHS,400*Nvar]);
set(gcf,'Name',[RunID ', ' varname{:} ' profiles']);
tiledlayout(Nvar, NPHS, 'TileSpacing', 'compact', 'Padding', 'compact');

for vi = 1:Nvar
    for fi = 1:NPHS
        nexttile;
        plot( squeeze(v{vi}(fi,:,opt.iPlt)), xp{vi});
        axis manual;
        hold on; plot(xlim, [0,0], 'k:', 'LineWidth', 1); hold off;
        
        ylim([min(xp{vi}), max(xp{vi})]);
        xlabel(varname{vi});
        
        % title shows phase names
        if vi==1, title([PHS{fi} '(t=0) = ' num2str(100*f0(fi), '%.0f'), '\%']); end
        
        % label z distance
        if fi==1
            ylabel(['$z$ [$\times ' zunits '$]']);
        else
            set(gca,'YTickLabel', []);
        end
        
        % legend?
        if vi==1 && fi==1 && opt.legd
            
            % make the legend showing times
            if max(t)<500, strspec = '%.0f';
            else, strspec = '%.1e';
            end
            leg = legend(num2str(t(opt.iPlt), strspec), 'location', 'best');
            title(leg, 't [s]');
        end
        
    end
end

if opt.save
    nvtmp = cell2mat(strcat(varname, '_'));
    SaveFigure([folder, RunID '/' RunID '_' nvtmp 'z_t' opt.fname], fig);
end


end

function [opt] = defopts (Nf,varargin)

% default opts
opt.fname  = '';        % extra filename info
opt.save   = 0;         % whether to save

opt.Nplt   = 1;         % number of lines to plot
opt.iPlt   = [];        % which timesteps to plot

opt.val    = 'mid';     % mid, mean, ind, what to plot
opt.xdsc   = 0;         % whether to divide x by max dsc
opt.xind   = 0;         % if val='ind', which x index to plot

opt.dfplt  = 0;         % if var == 'f', whether to plot change over time
opt.legd   = 1;         % legend showing time?

% allow structure alteration
args = reshape(varargin, 2, []);
for ia = 1:size(args,2)
    opt.(args{1,ia}) = args{2,ia};
end

% check form of Nplt
if isempty(opt.iPlt)
    if length(opt.Nplt)==1
        if opt.Nplt==1
            opt.iPlt = 1:Nf;
        else
            opt.iPlt = unique(round(linspace(1,Nf,opt.Nplt)));
        end
    else
        opt.iPlt   = opt.Nplt;
    end
end
opt.Nplt = length(opt.iPlt);

if opt.xind~=0
    opt.val = 'xind';
    opt.fname = ['_xi' num2str(opt.xind, '%.0f')];
end


end