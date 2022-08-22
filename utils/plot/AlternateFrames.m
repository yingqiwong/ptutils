function [fig,varmat,x] = AlternateFrames (folder, RunID, varname, varmat, varargin)
%
% fig = AlternateFrames (folder, RunID, varname, varmat, varargin)
% use this to make a looping gif between two frames of the same field
%
% EXAMPLES
% fig = AlternateFrames(folder, RunID, 'f');
%
%
% INPUTS
% folder    folder name where output folder is stored
% RunID     name of the run so that the total path to a mat file is
%               [folder  RunID / RunID_0.mat]
% varname   string of variable names
%           NB if you want certain axis limit options, varname has to be precise
% varmat    cell vector pair of variable matrices to plot
%               if mising or empty, load the variable of varname from file
%               varname, varmat must be same size
% varargin  plotting options (see defopts)
%
%
% OUTPUTS
% fig       figure handle
% varmat    matrix of variables plotted [NPHS x Nx x Nx]
% x         x positions
%
%
% DEFAULT OPTIONS
% opt.fname  = '';        % extra filename info
% opt.zzero  = 0;         % whether to plot a horizontal zero line
% opt.xdsc   = 0;         % whether to divide by initial max dsc
% opt.ti     = [];        % which time index
% opt.uaxes  = 1;         % whether to have uniform axes for all panels
% opt.Nstd   = 3;         % number of stds from mean for axis limits
%
% YQW, 14 June 2021
%

%  check inputs
Nvar = length(varname);
if nargin<4, varmat = []; end
if isempty(varmat) && ~iscell(varmat), varmat = cell(Nvar,1); end

% get plotting options
opt = defopts(varargin{:});

% load colormap
load('ocean.mat', 'ocean');

% get output mat files
[fp, fn] = GetOutputMatFiles(folder, RunID);
if isempty(opt.ti), opt.ti = [1,length(fn)]; end
load(fp, 'NPHS','D');

if opt.uaxes,  climits = zeros(NPHS, 2, Nvar); end

% load variables
[t, x, z, varmat] = LoadPlotVars(folder, RunID, varname, varmat);

if (opt.uaxes)
    climits = uniformaxislimits(opt.Nstd, varname, varmat, fp);
end


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
set(fig,'Name',RunID);
set(fig,'color','white');
hAx = default2dpanels(1, NPHS, 'aspectratio', length(x)/length(z), 'bot', 1.50);


% define filename
if ~isempty(opt.fname)
    fname = ['_', opt.fname];
else
    fname = varname;
end
filename = [folder, RunID '/' RunID '_altframe_' fname];


for j = 1:length(opt.ti)
    for iphs = 1:NPHS
        axes(hAx(iphs));
        
        imagesc(x,z,squeeze(varmat(iphs,:,:,opt.ti(j))));
        
        axis xy equal tight;
        if (opt.uaxes), caxis(climits(iphs,:)); end
        cb = colorbar; set(cb,TL{:},TS{:});
        
        hAx(iphs).YAxis.Exponent = 0;
        set(gca,TL{:},TS{:});
        if iphs>1
            set(gca,'YTickLabel',[]);
        else
            yl = ylabel(['depth [' zunit ']']);
            yl.Units = 'centimeters';
            yl.Position(1) = -1.2;
        end
        
        xlabel(['position [' zunit ']']);
        
        title(['$' varname '^',num2str(iphs),'$, t = ' num2str(t(opt.ti(j)),'%.1e') ' s'],TX{:},FS{:});
    end
    drawnow
    frame = getframe(fig);
    [A,map] = rgb2ind(frame2im(frame),256);
    
    if j == 1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.5);
    else
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.5);
    end
    
end

imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.5);



end



function [opt] = defopts (varargin)

% default opts
opt.fname  = '';        % extra filename info

opt.zzero  = 0;         % whether to plot a horizontal zero line
opt.xdsc   = 0;         % whether to divide by initial max dsc
opt.ti     = [];        % which time index

opt.uaxes  = 1;         % whether to have uniform axes for all panels
opt.Nstd   = 3;         % number of stds from mean for axis limits

% allow structure alteration
args = reshape(varargin, 2, []);
for ia = 1:size(args,2)
    opt.(args{1,ia}) = args{2,ia};
end


end









