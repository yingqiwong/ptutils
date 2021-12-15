function [fig,varmat,x] = PlayPhaseFracwTime (folder, RunID, varargin)
% 
% fig = PlayPhaseFracwTime (folder, RunID, varname, varmat, varargin)
% use this to make a videos showing the evolution of phase fraction in time. 
% accounted for box drift in periodic bcs
% Full 2D field.
%
% examples
% fig = PlayPhaseFracwTime(folder, RunID);
% fig = PlayPhaseFracwTime(folder, RunID, 'Nstd', 5);
% 
% 
% INPUTS
% folder    folder name where output folder is stored
% RunID     name of the run so that the total path to a mat file is
%               [folder  RunID / RunID_0.mat]
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
% YQW, 21 October 2021
%


% load colormap
load('../pantarhei/src/ocean.mat', 'ocean');

% get output mat files
fp = GetOutputMatFiles(folder, RunID);

load(fp, 'NPHS','BC','D');

if strcmp(BC, 'periodic')
    [t,x,z,f] = RmBoxDrift(folder, RunID);
else
    [t,x,f] = ExtractFieldwTime(folder, RunID, {'f'});
    x = x.*ones(length(t),1);
    z = x;
end

% get plotting options
opt = defopts(varargin{:});

if isempty(opt.iphs), opt.iphs = 1:NPHS; end
NPHS = length(opt.iphs);

if (opt.uaxes)
    climits = uniformaxislimits(opt.Nstd, 'f', f, fp);
end

Nf = length(t);

if (opt.xdsc)
    load(fp, 'delta0');
    x = x./max(delta0(:));
    z = z./max(delta0(:));
    zunit = 'dsc0';
else
    zunit = 'm';
    if floor(log10(D))>3
        % change units to km
        x = 1e-3*x; z = 1e-3*z; D = 1e-3*D;
        zunit = 'km';
    end
end

tyr = (365.25*24*60*60); tunit = 's';
if max(t)>tyr
    t = t/tyr; tunit = 'yr';
    if max(t)>1e3
        t = 1e-3*t; tunit = 'kyr';
    end
end


% prepare for plotting
TX = {'Interpreter','Latex'}; FS = {'FontSize',18};
TL = {'TickLabelInterpreter','Latex'}; TS = {'FontSize',14};


% initialize figure and axes
fig = figure;
colormap(ocean);
set(fig,'Position',[500,500,270*NPHS+50,250]);
set(fig,'Name',RunID);
set(fig,'color','white');

hAx = tight_subplot(1,NPHS,[0.05,0.04],[0.03,0.18],[0.08,0.04]);


% open movie file
filename = [folder, RunID '/' RunID '_fieldt_frmdrift'];
if ~isempty(opt.fname), filename = [filename '_' opt.fname]; end

vidObj = VideoWriter(filename, 'MPEG-4');
vidObj.FrameRate = opt.fps;
vidObj.Quality = 100;

open(vidObj);


for fi = 1:Nf
    for iplt = 1:NPHS
            axes(hAx(iplt));
            iphs = opt.iphs(iplt);
            
            imagesc(x(fi,:),z(fi,:),squeeze(f(iphs,:,:,fi)));
            shading interp
            
            axis xy equal tight;
            if (opt.uaxes), caxis(climits(iphs,:)); end
            cb = colorbar; set(cb,TL{:},TS{:});
            xlim([-D/2,D/2]);
            ylim([-D/2,D/2]);
            
            hAx(iplt).YAxis.Exponent = 0;
            set(gca,TL{:},TS{:}); set(gca,'XTickLabel',[]);
            if iplt>1
                set(gca,'YTickLabel',[]); 
            else
                yl = ylabel(['depth [' zunit ']']);
                yl.Position(1) = 1.2*min(x(1,:));
            end
            
            title(['$f^',num2str(iphs),'$'],TX{:},FS{:});
            
    end
    
    supertitle = ['t = ' num2str(t(fi),'%.2f') ' ' tunit];
    ha = annotation('textbox','Position',[0.5,0.99,0.05,0.02],...
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

opt.iphs   = [];
opt.zzero  = 0;         % whether to plot a horizontal zero line
opt.xdsc   = 0;         % whether to divide by initial max dsc

opt.uaxes  = 1;         % whether to have uniform axes for all panels
opt.Nstd   = 5;         % number of stds from mean for axis limits

opt.fps    = 15;         % frames per sec

% allow structure alteration
args = reshape(varargin, 2, []);
for ia = 1:size(args,2)
    opt.(args{1,ia}) = args{2,ia};
end


end









