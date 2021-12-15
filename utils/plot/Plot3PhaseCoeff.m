function fig = Plot3PhaseCoeff (f, K, varargin)
% 
% fig = Plot3PhaseCoeff (f, K, varargin)
% 
% example:  Plot3PhaseCoeff(f, cat(3,Kv,Kf), 'scl', 'log', 'PHS', PHS, 'cfname', {'K_v','K_f'});
% 
% plots the coefficients on a ternary diagram for a three phase system
% 
% INPUTS 
% f         phase fractions [NPHS x N]
% K         coefficients [NPHS x N x Nk]
% varargin  plotting and labeling options
% 
% YQW, 23 April 2021

% number of phases
NPHS = size(f,1);
np   = sqrt(size(f,2));

% number of coefficients
Ncf = size(K,3);

% options for plotting
opt = set_opts(NPHS, Ncf, varargin{:});

opt = getplotlimits(f, K, opt);


% assign cff and phs to row and columns
Nrow = Ncf;
Ncol = NPHS;


% load default plotting format
PlotDefFormat;

% prepare axes/borders dimensions
axh = 10;
axw = 8;
ahs = 2;
avs = 1.5;
axb = 0.2;
axt = 1.2;
axl = 1.4;
axr = 1.2;
fh = axb + Nrow*axh + (Nrow-1)*avs + axt;
fw = axl + Ncol*axw + (Ncol-1)*ahs + axr;


% initialize figure and axes
fig = figure;
set(fig,'defaultlinelinewidth',0.5);
set(fig,UN{[1,3]},'Position',[1 10 fw fh]);
set(fig,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
set(fig,'Color','w','InvertHardcopy','off');
set(fig,'Resize','off');

% set up axes
for irow = 1:Nrow
    for icol = 1:Ncol
        axInd = (irow-1)*Ncol + icol;
        ahpos = axl +    (icol-1)*axw +    (icol-1)*ahs;
        awpos = axb + (Nrow-irow)*axh + (Nrow-irow)*avs;
        ax(axInd) = axes(UN{[1,3]},'position',[ahpos awpos axw axh]);
    end
end



% prepare for plotting ternary plots
m = 5;
grids = linspace(0,1,m+1);
grids = grids(1:end-1);
labels = num2str(grids(2:end)');
[x3, y3] = terncoords(1-grids, grids, zeros(size(grids)));
[x2, y2] = terncoords(grids, zeros(size(grids)), 1-grids);
[x1, y1] = terncoords(zeros(size(grids)), 1-grids, grids);
n = m-1;
[X,Y]  =  terncoords(f(1,:),f(2,:),f(3,:));
tri    =  simpletri(np);

AxPhs = strcat('$\phi^{', opt.PHS(:), '}$');

% plot coefficients
for irow = 1:Nrow
    for icol = 1:Ncol
        axInd = (irow-1)*Ncol + icol;
        iphs = icol;
        icff = irow;
        
        axes(ax(axInd));
        
        % plot ternary axes
        hold on;
        set(gca,'visible','off');
        %maxz = max(log10(K(iphs,:,icff)));
        maxz = 1e10;
        plot3([0,1,0.5,0],[0,0,sin(1/3*pi),0],maxz.*ones(4,1),'color',[0,0,0],LW{[1,4]},'handlevisibility','off');
        axis equal tight;
        
        % plot grid on ternary axes
        for i = 1:n
            plot3([x1(i+1),x2(n-i+2)],[y1(i+1),y2(n-i+2)],maxz.*ones(2,1),'w','handlevisibility','off');
            plot3([x2(i+1),x3(n-i+2)],[y2(i+1),y3(n-i+2)],maxz.*ones(2,1),'w','handlevisibility','off');
            plot3([x3(i+1),x1(n-i+2)],[y3(i+1),y1(n-i+2)],maxz.*ones(2,1),'w','handlevisibility','off');
            plot3([x1(i+1),x2(n-i+2)],[y1(i+1),y2(n-i+2)],maxz.*ones(2,1),'k:','handlevisibility','off');
            plot3([x2(i+1),x3(n-i+2)],[y2(i+1),y3(n-i+2)],maxz.*ones(2,1),'k:','handlevisibility','off');
            plot3([x3(i+1),x1(n-i+2)],[y3(i+1),y1(n-i+2)],maxz.*ones(2,1),'k:','handlevisibility','off');
        end
        
        % plot coefficients
        switch opt.scl
            case 'linear'
                trisurf(tri,X,Y,K(iphs,:,icff));
                caxis(opt.cflim(:,icff));
            case 'log'
                trisurf(tri,X,Y,log10(K(iphs,:,icff)));
                caxis(log10(opt.cflim(:,icff)));
        end
        shading interp; view([0,90]);
        
        % label contours
        text(x3(2:end)+0.02,y3(2:end)-0.02,labels,FS{[1,3]},TX{:},HA{[1,2]},VA{[1,2]});
        text(x2(2:end)-0.00,y2(2:end)-0.02,labels,FS{[1,3]},TX{:},HA{[1,3]},VA{[1,4]});
        text(x1(2:end)-0.02,y1(2:end)-0.02,labels,FS{[1,3]},TX{:},HA{[1,4]},VA{[1,2]});
        
        % label phase names
        text( 1.02, 0.02,AxPhs{1}   ,TX{:},FS{[1,4]},UN{[1,2]},HA{[1,2]},VA{[1,3]});
        text( 0.05,-0.08,AxPhs{3}   ,TX{:},FS{[1,4]},UN{[1,2]},HA{[1,4]},VA{[1,3]});
        text(0.525,0.975,AxPhs{2}   ,TX{:},FS{[1,4]},UN{[1,2]},HA{[1,2]},VA{[1,2]});
        colormap(ax(axInd),ocean);
        
        % colorbar
        c = colorbar('location','southoutside',UN{[1,3]},FS{[1,3]},TL{:});
        cpos = get(c,'position');
        set(c,'position',cpos,TL{:});
        
        TitleText = ['\textbf{(' char(96+axInd) ')} ' opt.pltname{axInd}];
        text(-0.12,1.00,TitleText,TX{:},FS{[1,4]},UN{[1,2]},HA{[1,2]},VA{[1,3]});
        hold off; drawnow;
    end
end


end

function opt = set_opts (NPHS, Ncf, varargin)

opt.scl    = 'linear';
opt.PHS    = strcat({'f'}, num2str((1:NPHS)'));
opt.cfname = strcat({'K'}, num2str((1:Ncf)'));
opt.cflim  = [];

% allow structure alteration
args = reshape(varargin, 2, []);
for ia = 1:size(args,2)
    opt.(args{1,ia}) = args{2,ia};
end

% title of the plots
if ~isfield(opt, 'pltname')
    opt.pltname = cell(NPHS*Ncf,1);
    vi = 0;
    for ki = 1:Ncf
        for fi = 1:NPHS
            vi = vi+1;
            switch opt.scl
                case 'linear'
                    pre = '$';
                case 'log'
                    pre = 'log$_{10} ';
            end
            opt.pltname{vi} = strcat(pre, opt.cfname{ki}, '^{', opt.PHS{fi}, '}$');
        end
    end
end

end

function [opt] = getplotlimits (f, K, opt)

Ncf = size(K,3);

% collect limits of K so that we plot on uniform colorbar but not all
% extreme values
if isempty(opt.cflim)
    Npts = length(f(1,~isnan(f(1,:))));     % valid phase fractions
    knum = floor(0.05*Npts);                % 5th, 95th percentiles
    
    cflim = zeros(2,Ncf);
    for ki = 1:Ncf
        Kind = K(:,:,ki);
        
        % get extreme knum values
        minvals = mink(Kind(~isinf(Kind)&Kind>0), knum);
        maxvals = maxk(Kind(~isinf(Kind)&Kind>0), knum);
        
        cflim(1,ki) = minvals(end);
        cflim(2,ki) = maxvals(end);
    end
    opt.cflim = cflim;
end


end