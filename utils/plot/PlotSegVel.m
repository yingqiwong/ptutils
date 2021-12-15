function fig = PlotSegVel (PHS, f0Mat, fxi, w, wstar, wsegr)

NPHS = length(PHS);

fig(1) = figure;
set(gcf,'Position',[500,800,1000,350]);
lspec = {'LineStyle', 'none', 'MarkerSize', 12, 'linewidth', 3.5};
colororder(repmat(lines(2), 2, 1));

subplot(121);
for fi = 1:NPHS
        
    pos = w(fi,:)>=0;
    h = plot(f0Mat(fxi,pos), abs(w(fi,pos)), '+', lspec{:});
    uistack(h, 'bottom');
    
    hold on;
    plot(f0Mat(fxi,~pos), abs(w(fi,~pos)), '_', lspec{:});
end

set(gca,'yscale', 'log');

pos = wstar>=0;
semilogy(f0Mat(fxi,pos), abs(wstar(pos)), 'k+', lspec{:});
hstar = semilogy(f0Mat(fxi,~pos), abs(wstar(~pos)), 'k_', lspec{:});

hold off;

xlabel([PHS{fxi} ' fraction']);
ylabel('max $| w^i |$ [m/s]');
LegText = strcat('$w^{', PHS(:), '}$');
% legend([h, hstar], [LegText; {'$w^*$'}], 'Location', 'best');
title('(i) Vertical velocity');





subplot(122);
for fi = 1:NPHS    
    pos = wsegr(fi,:)>=0;
    h = plot(f0Mat(fxi,pos), abs(wsegr(fi,pos)), '+', lspec{:});
    uistack(h, 'bottom');

    hold on;
    plot(f0Mat(fxi,~pos), abs(wsegr(fi,~pos)), '_', lspec{:});
    
end
% hold off;
set(gca,'yscale', 'log');
hold off;

xlabel([PHS{fxi} ' fraction']);
ylabel('max $ |w^i_\Delta|$ [m/s]');
LegText = strcat('$w^{', PHS(:), '}_\Delta$');
legend(LegText, 'Location', 'best');
title('(ii) Vertical segregation velocity');





%{
fig(2) = figure;
[Nrow, Ncol] = GetSubplotRowCol(NPHS);
set(gcf,'Position',[500,200,500*Ncol,300*Nrow+50]);
for fi = 1:NPHS
    subplot(Nrow,Ncol,fi);
    plot(f0Mat(:,fxi), abs(wsegr(:,fi)./wstar), '+', lspec{:});
    xlabel([PHS{fxi} ' fraction']);
    ylabel('max $ |w^{i}_\Delta  /  w^*|  $');
    title(PHS{fi});
end
sgtitle('Vertical segr/reference velocity','FontSize',20);


if NPHS == 3
    addpath('../test/');

    fig(3) = figure;
    set(gcf,'Position',[500 500 1400 1100]);
    t = tiledlayout(3,4);
    
    
    for fi = 1:NPHS
        nexttile(fi);
        [X,Y] = PlotTernaryBG(f0Mat, PHS);
        scatter(X', Y', 100, w(:,fi), 'filled');
        hold on
        PlotColorbar(w(:,fi));
        text(-0.2,0.9,['$\langle w^{', PHS{fi}, '} \rangle$'],'FontSize',20);
        
        nexttile(fi+(NPHS+1));
        [X,Y] = PlotTernaryBG(f0Mat, PHS);
        scatter(X', Y', 100, wsegr(:,fi), 'filled');
        hold on
        PlotColorbar(wsegr(:,fi));
        text(-0.2,0.9,['$\langle w^{', PHS{fi}, '}_\Delta \rangle $'],'FontSize',20);
        
        nexttile(fi+2*(NPHS+1));
        [X,Y] = PlotTernaryBG(f0Mat, PHS);
        scatter(X', Y', 100, wsegr(:,fi)./abs(wstar), 'filled');
        hold on
        PlotColorbar(wsegr(:,fi)./abs(wstar));
        text(-0.2,0.9,['$\langle w^{', PHS{fi}, '}_\Delta \rangle / |\langle w^*\rangle| $'],'FontSize',20);
    end
    
    nexttile(NPHS+1);
    [X,Y] = PlotTernaryBG(f0Mat, PHS);
    scatter(X', Y', 100, wstar, 'filled');
    hold on
    PlotColorbar(wstar);
    text(-0.2,0.9,'$ \langle w^* \rangle$','FontSize',20);
    
    
    
end
%}

end

function [] = PlotColorbar (v)

v2   = log10(abs(v));
vOut = v(~isoutlier(v2) & ~isoutlier(v));
set(gca, 'CLim', [min(vOut), max(vOut)]);
c = colorbar('location','southoutside','TickLabelInterpreter','latex');
if all(vOut>=0) || all(vOut<=0)
    set(gca, 'ColorScale', 'log');
    ct = log10(abs(c.Limits));
    ct2 = floor(ct(1)):sign(c.Limits(1)):ceil(ct(2));
    c.Ticks = sort(sign(c.Limits(1))*10.^ct2,'ascend');
end
end


