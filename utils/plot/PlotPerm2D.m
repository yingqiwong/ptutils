function [figv, dsc, Kv, Kf, Cv, Cf, Xf] = PlotPerm2D (ffixi, ffix, eta0, d0, A, B, C, fxi, phsname)
% 
% [dsc, Kv, Kf, Cv, Cf, Xf] = PlotPerm2D (ffixi, ffix, eta0, d0, A, B, C)
% 
% Use this code to plot 2D lines of multi-phase permission functions,
% coefficients, seg-comp lengths. 
% 
% INPUTS
% ffixi     index of phase that you want to fix [NPHS-2 x 1]
%           NB: length of ffixi must be 2 less than total # of phases, 
%           this is just how the permutation goes
% ffix      values of phase that you want to fix [NPHS-2 x Nffix]
% eta0      pure-phase viscosity [NPHS x 1]
% d0        pure-phase grain scale [NPHS x 1]
% A,B,C     permission weight parameters for coefficient closure model
% fplti     index that you want to plot on the x axis [1x1] (optional)
% phsname   cell of phase names [NPHS x 1] (optional)

NPHS  = length(eta0);
Nffix = size(ffix,2);
Npts  = 401;

if nargin<8, fxi = 1; end
if nargin<9, phsname = num2str((1:NPHS)'); end


% check which phases are not fixed by input arguments
fvary = setdiff(1:NPHS, ffixi);

% connectivity figure
figv(1) = figure;
tiledlayout(NPHS,1,TileSpacing="compact",Padding="compact");
set(gcf,'defaultlinelinewidth',1);

% coefficient figure
figv(2) = figure;
tiledlayout(2,2,TileSpacing="compact",Padding="compact");
set(gcf,'defaultlinelinewidth',1);

% weights figure
figv(3) = figure;
tiledlayout(2,1,TileSpacing="compact",Padding="compact");
set(gcf,'defaultlinelinewidth',1);

figv(4) = figure;
set(gcf,'defaultlinelinewidth',1);
fpos = get(gcf,'Position'); set(gcf,'Position',[fpos(1),fpos(2)-500,fpos(3),fpos(4)+200]);

colors = lines(7);
lstyle = repmat({'-';'--';'-.';':'}, ceil(Nffix/4), 1);

for fi = 1:Nffix
    f1 = linspace(0, 1-sum(ffix(:,fi)), Npts);
    
    f = zeros(NPHS, Npts);
    f(ffixi,:) = ffix(:,fi).*ones(1,Npts);   % phases with fixed values
    f(fvary,:) = [f1; 1-f1-sum(ffix(:,fi))]; % phases that vary
    f(f<0) = nan;

    [dsc, Kv, Kf, Cv, Cf, Xf] = SegCompLength(f, eta0, d0, A, B, C);
    
    % plot connectivities
    figure(figv(1));
    for iphs = 1:NPHS
        nexttile(iphs);
        for jphs = 1:NPHS
            plot(f(fxi,:), squeeze(Xf(iphs,jphs,:)), lstyle{fi}, 'Color', colors(jphs,:)); 
            hold on;
        end
    end


    % plot coefficients
    figure(figv(2));
    nexttile(1); 
    for iphs=1:NPHS, semilogy(f(fxi,:), Kv(iphs,:)./f(iphs,:), lstyle{fi}, 'Color', colors(iphs,:)); hold on; end

    nexttile(2); 
    for iphs=1:NPHS, semilogy(f(fxi,:), Kf(iphs,:)./f(iphs,:), lstyle{fi}, 'Color', colors(iphs,:)); hold on; end
    
    nexttile(3);     
    for iphs=1:NPHS, semilogy(f(fxi,:), f(iphs,:).^2./Cv(iphs,:), lstyle{fi}, 'Color', colors(iphs,:)); hold on; end
    
    nexttile(4);
    for iphs=1:NPHS, semilogy(f(fxi,:), f(iphs,:).^2./Cf(iphs,:), lstyle{fi}, 'Color', colors(iphs,:)); hold on; end

    
    % plot pressure weights, velocity weights
    figure(figv(3));
    nexttile(1);
    for iphs=1:NPHS, plot(f(fxi,:), Cf(iphs,:)./sum(Cf), lstyle{fi}, 'Color', colors(iphs,:)); hold on; end

    nexttile(2);
    for iphs=1:NPHS, plot(f(fxi,:), Cv(iphs,:)./sum(Cv), lstyle{fi}, 'Color', colors(iphs,:)); hold on; end
    
    % plot dsc
    figure(figv(4));
    semilogy(f(fxi,:), squeeze(dsc(1,2,:)), lstyle{fi}, 'Color', colors(1,:)); hold on; 
    semilogy(f(fxi,:), squeeze(dsc(1,3,:)), lstyle{fi}, 'Color', colors(2,:)); hold on; 
    semilogy(f(fxi,:), squeeze(dsc(2,3,:)), lstyle{fi}, 'Color', colors(3,:)); hold on; 
    semilogy(f(fxi,:), squeeze(dsc(3,2,:)), lstyle{fi}, 'Color', colors(4,:)); hold on; 


end


% finish up figure formatting
figure(figv(1));
for iphs = 1:NPHS
    nexttile(iphs); ylabel('Connectivity'); title(['$\phi^{' phsname{iphs} '}$']);
end
legend(strcat('$', phsname, '$'), Location='east'); 
xlabel(['$\phi^{' phsname{fxi} '}$']);

legtext = strcat('$', repmat(phsname(:),Nffix,1), {'$, '}, repelem(num2str(ffix(:)),NPHS,1));

figure(figv(2));
nexttile(1); title('$K_v/\phi$'); 
nexttile(2); title('$K_f/\phi$');
nexttile(3); title('$\phi^2/C_v$'); xlabel(['$\phi^{' phsname{fxi} '}$']);
nexttile(4); title('$\phi^2/C_f$'); xlabel(['$\phi^{' phsname{fxi} '}$']);
legend(legtext, Location='best',NumColumns=Nffix, Interpreter='latex'); 


figure(figv(3));
nexttile(1); title('Pressure weights'); 
legend(legtext, Location='best',NumColumns=Nffix, Interpreter='latex'); 
nexttile(2); title('Velocity weights');
xlabel(['$\phi^{' phsname{fxi} '}$']);

figure(figv(4));
hold on; plot(xlim, d0(1)*ones(1,2), 'k:'); hold off;
ylimits = ylim; ylim([1e-10,ylimits(2)]);
title('seg-comp length [m]') ; xlabel(['$\phi^{' phsname{fxi} '}$']);
legtext = strcat(repmat(...
                [strcat(phsname{2},'-',phsname{1},', '); ...
                 strcat(phsname{3},'-',phsname{1},', '); ...
                 strcat(phsname{3},'-',phsname{2},', '); ...
                 strcat(phsname{2},'-',phsname{3},', ')],NPHS,1), ...
             repelem(num2str(ffix(:)),4,1));
legend(legtext, Location='southoutside',NumColumns=Nffix);


end