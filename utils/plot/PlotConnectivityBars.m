function [fig] = PlotConnectivityBars (f, Xf, varargin)
% 
% PlotConnectivityBars (f, Xf, varargin)
% 
% example:  PlotConnectivityBars (f, Xf, 'folder', '../out/')
% 
% INPUTS
% f         phase fractions [NPHS x Nplt]
% Xf        permission weights [NPHS x NPHS x Nplt]
% varargin  plotting options
% 
% YQW, 6 May 2021

% get lengths
[NPHS, Nplt] = size(f);

% get plotting options
opt  = defopts(NPHS,varargin{:});

fig=figure;
set(gcf,'Position',[500,500,400,100*NPHS]);
colors = lines(NPHS);

i = 0; inds = []; 
for xfi = 1:Nplt
    
    % plot the stacked bars of connectivities
    for iphs = 1:NPHS
        h1 = barh(i+1, Xf(iphs,:,xfi), 'stacked'); hold on;
        
        % assign colors to the bars
        for iphs2 = 1:NPHS
            h1(iphs2).FaceColor = colors(iphs2,:);
        end
        
        text(1.05,i+1,num2str(100*f(iphs,xfi),'%.0f'),'FontSize',16);
        inds = [inds, i+1];

        i = i+1;
    end
    i = i+3;
end

xlim([0,1.15]);
barlabs = repmat(opt.phasenames(:), size(Xf,3), 1);
set(gca,'YTick',inds,'YTickLabel',barlabs,'YDir','reverse');
title('Connectivity $X_\phi$');

if opt.save
    SaveFigure([opt.folder, 'connectivity', opt.fname]);
end


end

function [opt] = defopts (NPHS,varargin)
% default opts

opt.phasenames = strcat({'f'}, num2str((1:NPHS)', '%d'));

% saving options
opt.save   = 1;        % whether to save
opt.folder = '';       % folder to save figure to
opt.fname  = '';       % extra filename info


% allow structure alteration
args = reshape(varargin, 2, []);
for ia = 1:size(args,2)
    opt.(args{1,ia}) = args{2,ia};
end

end