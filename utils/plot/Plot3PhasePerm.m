function ax = Plot3PhasePerm (f, Xf, PHS)
% 
% fig = Plot3PhaseCoeff (f, K, varargin)
% 
% example:  Plot3PhaseCoeff(f, cat(3,Kv,Kf), 'scl', 'log', 'PHS', PHS, 'cfname', {'K_v','K_f'});
% 
% plot connectivities of three phase systems
% 
% INPUTS 
% f         phase fractions [NPHS x N]
% Xf        coefficients [NPHS x NPHS x N x Nk]
% varargin  plotting and labeling options
% 
% YQW, 23 April 2021


% set up
NPHS = size(f,1);
np   = sqrt(size(f,2));
if nargin<3, PHS = strcat({'f'}, num2str((1:NPHS)')); end

caxmax = inf;
if max(Xf,[],'all') == 1, caxmax = 1; end

PlotDefFormat;

% prepare axes/borders dimensions
axh = 8;
axw = 8;
ahs = 2.5;
avs = 1.5;
axb = 1.8;
axt = 1.2;
axl = 1.4;
axr = 1.2;
fh = axb + NPHS*axh + (NPHS-1)*avs + axt;
fw = axl + NPHS*axw + (NPHS-1)*ahs + axr;


% initialize figure and axes
fig = figure;
set(fig,'defaultlinelinewidth',0.5);
set(fig,UN{[1,3]},'Position',[1 10 fw fh]);
set(fig,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
set(fig,'Color','w','InvertHardcopy','off');
set(fig,'Resize','off');

% set up axes
ax(1) = axes(UN{[1,3]},'position',[axl+0*axw+0*ahs axb+2*axh+2*avs axw axh]);
ax(2) = axes(UN{[1,3]},'position',[axl+1*axw+1*ahs axb+2*axh+2*avs axw axh]);
ax(3) = axes(UN{[1,3]},'position',[axl+2*axw+2*ahs axb+2*axh+2*avs axw axh]);
ax(4) = axes(UN{[1,3]},'position',[axl+0*axw+0*ahs axb+1*axh+1*avs axw axh]);
ax(5) = axes(UN{[1,3]},'position',[axl+1*axw+1*ahs axb+1*axh+1*avs axw axh]);
ax(6) = axes(UN{[1,3]},'position',[axl+2*axw+2*ahs axb+1*axh+1*avs axw axh]);
ax(7) = axes(UN{[1,3]},'position',[axl+0*axw+0*ahs axb+0*axh+0*avs axw axh]);
ax(8) = axes(UN{[1,3]},'position',[axl+1*axw+1*ahs axb+0*axh+0*avs axw axh]);
ax(9) = axes(UN{[1,3]},'position',[axl+2*axw+2*ahs axb+0*axh+0*avs axw axh]);



% prepare for plotting ternary plots
[X,Y]  =  terncoords(f(1,:),f(2,:),f(3,:));
tri    =  simpletri(np);

AxPhs = strcat('$\phi^{', PHS(:), '}$');

axInd = 0;
for k = 1:NPHS
    for j = 1:NPHS
        axInd = axInd + 1;
        
        axes(ax(axInd));
        terngrid(5);

        % plot connectivities
        trisurf(tri,X,Y,squeeze(Xf(k,j,:)));
        shading interp; view([0,90]);
        vertexlabel(AxPhs{1},AxPhs{2},AxPhs{3});
        
        % colorbar
        colormap(ax(axInd),ocean);

        if caxmax == 1
            caxis([0,caxmax]);
            if axInd == 9
                c = colorbar('location','southoutside',UN{[1,3]},FS{[1,3]},TL{:},'Ticks',[0:0.2:1]);
                cpos = get(c,'position'); cpos(1) = axl+2*axw+2*ahs; cpos(2) = axb-0.50*avs; cpos(3) = axw; cpos(4) = 0.4;
                set(c,'position',cpos,TL{:});
            end
        else
            c = colorbar('location','southoutside',UN{[1,3]},FS{[1,3]},TL{:});
        end
        
            
        TitleText = ['\textbf{(' char(96+axInd) ')}~$X_\phi^{' PHS{j} '-' PHS{k} '}$' ];
        text(-0.12,1.00,TitleText,TX{:},FS{[1,4]},UN{[1,2]},HA{[1,2]},VA{[1,3]});
        hold off; drawnow;
    end
end
end
