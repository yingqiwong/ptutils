function [fig,varmat,x] = PlayFieldwTime (folder, RunID, varname, varmat, varargin)
% 
% fig = PlotFieldwTime (folder, RunID, varname, varmat, varargin)
% use this to make a videos showing the evolution of a field in time. 
% Full 2D field.
%
% examples
% fig = PlotFieldwTime(folder, RunID, 'f');
% fig = PlotFieldwTime(folder, RunID, 'f', [], 'Nstd', 5);
% fig = PlotFieldwTime(folder, RunID, 'rhobar', rhobar, 'Nstd', 5);
% 
% 
% INPUTS
% folder    folder name where output folder is stored
% RunID     name of the run so that the total path to a mat file is
%               [folder  RunID / RunID_0.mat]
% varname   cell vector of variable names
%           NB if you want certain axis limit options, varname has to be precise
% varmat    cell vector of variable matrices to plot
%               if mising or empty, load the variable of varname from file
%               varname, varmat must be same size
% varargin  plotting options (see defopts)
%
% OUTPUTS
% fig       figure handle
% varmat    matrix of variables plotted [NPHS x Nx x Nx]
% x         x positions
%
% DEFAULT OPTIONS
% opt.fname  = '';        % extra filename info
% opt.zzero  = 0;         % whether to plot a horizontal zero line
% opt.xdsc   = 0;         % whether to divide by initial max dsc
% opt.uaxes  = 1;         % whether to have uniform axes for all panels
% opt.Nstd   = 3;         % number of stds from mean for axis limits
% 
% YQW, editted 4 May 2021
%

%  check inputs
if nargin<4, varmat = cell(length(varname),1); end

% load colormap
load('ocean.mat', 'ocean');

% get output mat files
fp = GetOutputMatFiles(folder, RunID);

load(fp, 'NPHS','D');
Nvar = length(varname);

% get plotting options
opt = defopts(varargin{:});

if opt.uaxes, climits = zeros(NPHS, 2, Nvar); cblimits = climits; end

% load variables
for vi = 1:Nvar
    [t, x, z, varmat{vi}] = LoadPlotVars(folder, RunID, varname{vi}, varmat{vi});
    
    if (opt.uaxes)
        [climits(:,:,vi),cblimits(:,:,vi)] = uniformaxislimits(opt.Nstd, varname{vi}, varmat{vi}, fp);
    end
end

Nf = size(varmat{1},4);

if (opt.xdsc)
    load(fp, 'delta0');
    x = x./max(delta0(:));
    z = z./max(delta0(:));
    zunit = 'dsc0';
else
    zunit = 'm';
    if floor(log10(D(1)))>3
        % change units to km
        x = 1e-3*x; z = 1e-3*z; D = 1e-3*D;
        zunit = 'km';
    end
end


% prepare for plotting
TX = {'Interpreter','Latex'}; FS = {'FontSize',18};
TL = {'TickLabelInterpreter','Latex'}; TS = {'FontSize',14};


% initialize figure and axes
fig = figure;
colormap(ocean);
set(fig,'Position',[500,500,300*NPHS+50,240*Nvar+60]);
set(fig,'Name',RunID);
set(fig,'color','white');
hAx = default2dpanels(1, NPHS, 'aspectratio', length(x)/length(z), 'bot', 1.50);

% open movie file
if ~isempty(opt.fname)
    fname = ['_', opt.fname]; 
else
    fname = strjoin(varname,'');
end
filename = [folder, RunID '/' RunID '_fieldt_' fname];

vidObj = VideoWriter(filename, 'MPEG-4');
vidObj.FrameRate = opt.fps;
vidObj.Quality = 100;

open(vidObj);


for fi = 1:Nf
    for iphs = 1:NPHS
        for ivar = 1:Nvar
            axes(hAx((ivar-1)*NPHS+iphs));
            
            imagesc(x,z,squeeze(varmat{ivar}(iphs,:,:,fi)));
            
            axis xy equal tight;
            if (opt.uaxes), caxis(climits(iphs,:,ivar)); end
            cb = colorbar; set(cb,TL{:},TS{:});
            cb.Limits = cblimits(iphs,:,ivar);
            
            if iphs>1
                set(gca,'YTickLabel',[]);
            else
                yl = ylabel(['depth [' zunit ']']);
                yl.Units = 'centimeters';
                yl.Position(1) = -1.2;
            end
            
            xlabel(['position [' zunit ']']);
            
            hAx((iphs-1)*Nvar+ivar).YAxis.Exponent = 0;
            set(gca,TL{:},TS{:}); 
            if iphs>1, set(gca,'YTickLabel',[]); end
            
            title(['$' varname{ivar} '^',num2str(iphs),'$'],TX{:},FS{:});
            
        end
    end
    
    supertitle = ['t = ' num2str(t(fi),'%.1e') ' s'];
    ha = annotation('textbox','Position',[0.5,0.98,0.05,0.02],...
        'String', supertitle, ...
        'HorizontalAlignment','center','VerticalAlignment','top',...
        'FitBoxToText','on','EdgeColor','none','FontSize',22,TX{:});
    
    currFrame = getframe(fig);
    writeVideo(vidObj,currFrame);
    
    ha.String = '';
end

close(vidObj);

end



function [opt] = defopts (varargin)

% default opts
opt.fname  = '';        % extra filename info

opt.zzero  = 0;         % whether to plot a horizontal zero line
opt.xdsc   = 0;         % whether to divide by initial max dsc

opt.uaxes  = 1;         % whether to have uniform axes for all panels
opt.Nstd   = 3;         % number of stds from mean for axis limits

opt.fps    = 20;         % frames per sec

% allow structure alteration
args = reshape(varargin, 2, []);
for ia = 1:size(args,2)
    opt.(args{1,ia}) = args{2,ia};
end


end









