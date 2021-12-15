function [s, e1, e2] = CalcPrincipalStress_old (folder, RunID, ti, plot_opt, bgname, bgmat, t, x, varargin)
%
% [s, e1, e2] = CalcPrincipalStress (folder, RunID)
%
% calculates the principal stresses and principal stress directions at each
% grid point (i.e. eigenvalues and eigenvectors of stress matrix).
% 2 principal stresses and 2 directions because 2D.
% easy to calculate analytically.
%
% 
% INPUTS
% folder    folder name where output folder is stored
% RunID     name of the run so that the total path to a mat file is
%               [folder  RunID / RunID_0.mat]
%
% 
% OUTPUTS
% s         matrix of principal stresses [NPHS x N x N x Nt x 2]
% e1        1st stress direction (corresponding to s1) [NPHS x N x N x Nt x 2]
% e2        2nd stress direction (corresponding to s2) [NPHS x N x N x Nt x 2]
%
% 
% DEFAULT OPTIONS
% opt.fname  = '';        % specified filename
% opt.save   = 0;         % whether to save
% 
% opt.Nquiv  = 21;        % number of quiver arrows (downsampling)
% opt.bcind  = 0.05;      % exclude arrows from boundary
% opt.uquiv  = 1;         % whether to plot same length quiver for all plts
% opt.splt   = 1;         % whether to plot 1 or 2 principal stresses
% 
% opt.zline  = [];        % whether to plot a horizontal zero line
% opt.xdsc   = 0;         % whether to divide by initial max dsc
% opt.boxind = [];        % whether to plot a subset of domain
% 
% opt.Nplt   = 5;         % number of panels to plot
% opt.iPlt   = [];        % which panels to plot
% opt.uaxes  = 1;         % whether to have uniform axes for all panels
% opt.Nstd   = 3;         % number of stds from mean for axis limits
% 
% YQW, 8 June 2021

if isempty(ti)
    % load components of stress tensor
    [t, x, varmat] = ExtractFieldwTime(folder, RunID, {'qvxx','qvxz','qvzz','f','p'});
else
    [~,fn] = GetOutputMatFiles(folder, RunID);
    varmat = struct2cell(load(fn{ti}, 'qvxx','qvxz','qvzz','f','p'));
    varargin = [varargin, 'iPlt', 1];
    
    if nargin<7, load(fn{ti}, 'x','time'); t = time; end
end

% take negative of qv to get positive tension stress matrix
% remove f x p to get deviatoric component only
qvxx = - (varmat{1} - varmat{4}.*varmat{5});
qvxz = -  varmat{2};
qvzz = - (varmat{3} - varmat{4}.*varmat{5});

% interpolate shear stresses onto pressure grid (cell centers)
qvxzcc = (qvxz(:,1:end-1,1:end-1,:) + qvxz(:,1:end-1,2:end,:) + qvxz(:,2:end,1:end-1,:) + qvxz(:,2:end,2:end,:)) / 4;

% calculate principal stresses
s1 = 0.5*(qvxx+qvzz) + sqrt(0.25*(qvxx-qvzz).^2 + qvxzcc.^2);
s2 = 0.5*(qvxx+qvzz) - sqrt(0.25*(qvxx-qvzz).^2 + qvxzcc.^2);
s  = cat(5, s1, s2);

% calculate first eigenvector (normalised)
e11 = -qvxzcc./(qvxx-s1);
e1  = 1./sqrt(1+e11.^2).*cat(5, e11, ones(size(e11)));

% calculate second eigenvector (normalised)
e21 = -qvxzcc./(qvxx-s2);
e2  = 1./sqrt(1+e21.^2).*cat(5, e21, ones(size(e21)));

if nargin < 4, return; end
if nargin < 5, bgname = 'fp'; end
if nargin < 6, bgmat = []; end
if isempty(bgname), bgname = 'fp'; end
switch bgname
    case 's1'   , bgmat = s1; 
    case 's2'   , bgmat = s2; 
    case 'fp'   , bgmat = varmat{4}.*varmat{5};
    case 'sdiff', bgmat = s1-s2;
end
if plot_opt, plot_pstress(folder, RunID, bgname, bgmat, t, x, s, e1, e2, varargin{:}); end

end



function plot_pstress (folder, RunID, bgname, bg, t, x, s, e1, e2, varargin)

fp = GetOutputMatFiles(folder, RunID);
Nf = size(s,4);

% get plotting options
opt = defopts(Nf, varargin{:});

% load colormap
load('../pantarhei/src/ocean.mat', 'ocean');

% load plotting variables
[t, x, bg] = LoadPlotVars(folder, RunID, {bgname}, bg, t, x);

smax = max(abs(s),[],[5]);
se1 = (1+s(:,:,:,:,1)./smax).*e1;
se2 = (1+s(:,:,:,:,2)./smax).*e2;

if (opt.xdsc)
    load(fp, 'delta0');
    x  = x./max(delta0(:));
end

% check sizes and make sure bg, s, e are the same size
NPHS = size(se1,1);
if size(bg,1)<NPHS, bg = bg.*ones(NPHS,1); end

% check if want to plot a subset of domain
[xplt, zplt, bg, se1, se2] = domainsubset(opt.boxind, x, bg, se1, se2);

% prepare quiver
scl = 4*median(abs([se1,se2]),[2,3,4,5],'omitnan');
[Xquiv, Zquiv, se1, se2] = quivds(xplt, zplt, opt, scl, se1, se2);


% get axis limits
if opt.uaxes, climits = uniformaxislimits(opt.Nstd, bgname, bg, fp); end

% prepare for plotting
TX = {'Interpreter','Latex'}; FS = {'FontSize',18};
TL = {'TickLabelInterpreter','Latex'}; TS = {'FontSize',14};

% initialize figure and axes
fig = figure;
colormap(ocean);
set(fig,'Name',[RunID '_stressfield']);

axh = 240; axb =  3; axgh =  3; axt = 40;
axw = 260; axl = 50; axgw = 10; axr = 20; 
fh = axb +     NPHS*axh + (    NPHS-1)*axgh + axt;
fw = axl + opt.Nplt*axw + (opt.Nplt-1)*axgw + axr;
set(fig,'Position',[500,500,fw,fh]);
hAx = tight_subplot(NPHS,opt.Nplt,[axgh/fh,axgw/fw], [axb,axt]/fh, [axl,axr]/fw);

% plot panels
for m = 1:opt.Nplt
    
    for iphs=1:NPHS
        axes(hAx((iphs-1)*opt.Nplt+m));
        
        % plot background
        imagesc(xplt,zplt,squeeze(bg(iphs,:,:,opt.iPlt(m))));
        hold on;
        
        % plot zero line
        if ~isempty(opt.zline)        % plot a horizontal line
            zloc = x(opt.zline);
            plot(xlim'.*ones(1,length(opt.zline)), zloc.*ones(2,1), 'LineWidth',0.5,'Color',0.5*ones(1,3));
        end
        
        if any(opt.splt == 2)
            hquiv = quiver(Xquiv,Zquiv,squeeze(se2(iphs,:,:,opt.iPlt(m),1)),  squeeze(se2(iphs,:,:,opt.iPlt(m),2)),'color', 0.5*ones(1,3),'LineWidth',1);
            hquiv.AutoScale = 'off'; hquiv.ShowArrowHead = 'off';
            hquiv = quiver(Xquiv,Zquiv,-squeeze(se2(iphs,:,:,opt.iPlt(m),1)),-squeeze(se2(iphs,:,:,opt.iPlt(m),2)),'color', 0.5*ones(1,3),'LineWidth',1);
            hquiv.AutoScale = 'off'; hquiv.ShowArrowHead = 'off';
        end
        
        if any(opt.splt==1)
            hquiv = quiver(Xquiv,Zquiv,squeeze(se1(iphs,:,:,opt.iPlt(m),1)),squeeze(se1(iphs,:,:,opt.iPlt(m),2)),'k','LineWidth',1);
            hquiv.AutoScale = 'off'; hquiv.ShowArrowHead = 'off';
            hquiv = quiver(Xquiv,Zquiv,-squeeze(se1(iphs,:,:,opt.iPlt(m),1)),-squeeze(se1(iphs,:,:,opt.iPlt(m),2)),'k','LineWidth',1);
            hquiv.AutoScale = 'off'; hquiv.ShowArrowHead = 'off';
        end

        hold off;
        
        axis xy equal tight;
        xlim([min(xplt), max(xplt)]); ylim([min(zplt), max(zplt)]);
        if (opt.uaxes), caxis(climits(iphs, :)); end
        
        cb = colorbar; set(cb,TL{:},TS{:});
        hAx((iphs-1)*opt.Nplt+m).YAxis.Exponent = 0;
        set(gca,TL{:},TS{:}); set(gca,'XTickLabel',[]);
        if m>1, set(gca,'YTickLabel',[]); end
        
        if iphs==1, title(['t = ' num2str(t(opt.iPlt(m)),'%.1e') ' s'],TX{:},FS{:}); end
        
    end
end



% plot supertitle to tell us what is plotted
supertitle = ['bg: $' bgname '$; vectors: principal stresses'];
ha = annotation('textbox','Position',[0.5,0.98,0.05,0.02],'String', supertitle, ...
    'HorizontalAlignment','center','VerticalAlignment','top',...
    'FitBoxToText','on','EdgeColor','none','FontSize',22,TX{:});



if opt.save
    if ~isempty(opt.fname)
        fname = opt.fname;
    else
        fname = ['s' num2str(1:opt.splt, '%.0f')];
    end
    figname = [folder, RunID '/' RunID '_pstresst_' fname];
    SaveFigure(figname, fig);
end


end







function [opt] = defopts (Nf, varargin)

% default opts
opt.fname  = '';        % specified filename
opt.save   = 0;         % whether to save

opt.Nquiv  = 15;        % number of quiver arrows (downsampling)
opt.bcind  = 0.05;      % exclude arrows from boundary
opt.uquiv  = 0;         % whether to plot same length quiver for all plts
opt.splt   = [1,2];         % whether to plot 1 or 2 principal stresses

opt.zline  = [];        % whether to plot a horizontal zero line
opt.xdsc   = 0;         % whether to divide by initial max dsc
opt.boxind = [];        % whether to plot a subset of domain

opt.Nplt   = 5;         % number of panels to plot
opt.iPlt   = [];        % which panels to plot
opt.uaxes  = 1;         % whether to have uniform axes for all panels
opt.Nstd   = 3;         % number of stds from mean for axis limits


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






