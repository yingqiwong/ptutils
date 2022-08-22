function [fig,varmat] = PlotVarDistribwTime (folder, RunID, varname, varmat, varargin)
%
% fig = PlotVarDistribwTime (folder, RunID, varname, varmat, varargin)
% use this to plot the distribution of variables through time
%
% examples
% fig = PlotVarDistribwTime(folder, RunID, 'f');
% fig = PlotVarDistribwTime(folder, RunID, 'rhobar', rhobar, 'uaxes', 0);
%
%
% INPUTS
% folder    folder name where output folder is stored
% RunID     name of the run so that the total path to a mat file is
%               [folder  RunID / RunID_0.mat]
% varname   character vector of the name of variable to plot
%           NB if you want certain axis limit options, varname has to be precise
% varmat    matrix of the variable to plot
%               if mising or empty, load the variable of varname from file
% varargin  plotting options (see defopts)
%
% OUTPUTS
% fig       figure handle
% varmat    matrix of variables plotted [NPHS x Nx x Nx]
% x         x positions
%
% YQW, editted 4 May 2021
%

%  check inputs
if nargin<4, varmat = []; end

% collect varmat to plot
[t, ~, ~, varmat] = LoadPlotVars(folder, RunID, varname, varmat);
Nf = length(t);

% get plotting options
opt  = defopts(Nf, varargin{:});
NPHS = size(varmat, 1);


% get uniform axes limits?
if (opt.uaxes)
    % limits = range of values
    climits = [min(varmat, [], [2,3,4]), max(varmat, [], [2,3,4])];
end

% get histogram bin edges
edges = zeros(NPHS, opt.Nbins);
for iphs = 1:NPHS
    edges(iphs,:) = linspace(climits(iphs,1), climits(iphs,2), opt.Nbins);
end



% prepare for plotting
TX = {'Interpreter','Latex'}; FS = {'FontSize',18};
TL = {'TickLabelInterpreter','Latex'}; TS = {'FontSize',14};


% initialize figure and axes
fig = figure;
set(fig,'Name',RunID);

axh = 240; axb = 30; axgh =  3; axt = 30;
axw = 260; axl = 50; axgw = 40; axr = 20;
fh = axb +     NPHS*axh + (    NPHS-1)*axgh + axt;
fw = axl + opt.Nplt*axw + (opt.Nplt-1)*axgw + axr;
set(fig,'Position',[500,500,fw,fh]);
hAx = tight_subplot(NPHS,opt.Nplt,[axgh/fh,axgw/fw], [axb,axt]/fh, [axl,axr]/fw);

ylimits = zeros(NPHS, opt.Nplt, 2);

% plot panels
for m = 1:opt.Nplt
    
    for iphs = 1:NPHS
        axes(hAx((iphs-1)*opt.Nplt+m));
        
        if m>1
            histogram(squeeze(varmat(iphs,:,:,opt.iPlt(1))), edges(iphs,:), 'EdgeColor', 'none');
            hold on;
        end
        
        histogram(squeeze(varmat(iphs,:,:,opt.iPlt(m))), edges(iphs,:), 'EdgeColor', 'none');
%         set(gca,'YScale','log');
        if (opt.uaxes), xlim(climits(iphs, :)); ylimits(iphs,m,:) = ylim; end
        title(['$' varname '^',num2str(iphs),'$, t = ' num2str(t(opt.iPlt(m)),'%.1e') ' s'],TX{:},FS{:});
    end
end

% now go through and make all y axes uniform
if (opt.uaxes)
    for iphs = 1:NPHS
        ylimfinal = [1, max(ylimits(iphs,:,2),[],2)];
        for m = 1:opt.Nplt
            axes(hAx((iphs-1)*opt.Nplt+m)); ylim(ylimfinal(1,:));
        end
    end
end

if opt.save
    if ~isempty(opt.fname), opt.fname = ['_', opt.fname]; end
    SaveFigure([folder, RunID '/' RunID '_varhist_' varname opt.fname], fig);
end



end



function [opt] = defopts (Nf, varargin)

% default opts
opt.fname  = '';       % extra filename info
opt.save   = 0;        % whether to save
opt.Nplt   = 3;        % number of panels to plot
opt.iPlt   = [];

opt.Nbins  = 100;      % number of histogram bins
opt.uaxes  = 1;        % whether to have uniform axes for all panels

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


end