function [dsc, Kv, Kf, Cv, Cf, Xf] = PlotPerm2D (ffixi, ffix, eta0, d0, A, B, C, fplti)
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

if nargin<8, fplti = 1; end

NPHS  = length(eta0);
Nffix = size(ffix,2);
Npts  = 401;

% check which phases are not fixed by input arguments
fvary = setdiff(1:NPHS, ffixi);

% connectivity figure
figXf = figure;
tiledlayout(NPHS,1,TileSpacing="compact",Padding="compact");
set(gcf,'defaultlinelinewidth',1);

% coefficient figure
figKC = figure;
tiledlayout(2,2,TileSpacing="compact",Padding="compact");
set(gcf,'defaultlinelinewidth',1);

% weights figure
figWgt = figure;
tiledlayout(NPHS,1,TileSpacing="compact",Padding="compact");
set(gcf,'defaultlinelinewidth',1);

colors = lines(NPHS);


for fi = 1:Nffix
    f1 = linspace(0, 1-sum(ffix(:,fi)), Npts);
    
    f = zeros(NPHS, Npts);
    f(ffixi,:) = ffix(:,fi).*ones(1,Npts);   % phases with fixed values
    f(fvary,:) = [f1; 1-f1-sum(ffix(:,fi))]; % phases that vary
    f(f<0) = nan;

    [dsc, Kv, Kf, Cv, Cf, Xf] = SegCompLength(f, eta0, d0, A, B, C);
    
    % plot connectivities
    figure(figXf);
    for iphs = 1:NPHS
        nexttile(iphs);
        for jphs = 1:NPHS
            plot(f(fplti,:), squeeze(Xf(iphs,jphs,:)), 'Color', colors(jphs,:).^(1/fi)); 
            hold on;
        end
    end


    % plot coefficients
    figure(figKC);
    nexttile(1); 
    for iphs=1:NPHS, semilogy(f(fplti,:), Kv(iphs,:)./f(iphs,:), 'Color', colors(iphs,:).^(1/fi)); hold on; end

    nexttile(2); 
    for iphs=1:NPHS, semilogy(f(fplti,:), Kf(iphs,:)./f(iphs,:), 'Color', colors(iphs,:).^(1/fi)); hold on; end
    
    nexttile(3);     
    for iphs=1:NPHS, semilogy(f(fplti,:), f(iphs,:).^2./Cv(iphs,:), 'Color', colors(iphs,:).^(1/fi)); hold on; end
    
    nexttile(4);
    for iphs=1:NPHS, semilogy(f(fplti,:), f(iphs,:).^2./Cf(iphs,:), 'Color', colors(iphs,:).^(1/fi)); hold on; end

    
    % plot pressure weights, velocity weights, dsc
    figure(figWgt);
    nexttile(1);
    for iphs=1:NPHS, plot(f(fplti,:), Cf(iphs,:)./sum(Cf), 'Color', colors(iphs,:).^(1/fi)); hold on; end

    nexttile(2);
    for iphs=1:NPHS, plot(f(fplti,:), Cv(iphs,:)./sum(Cv), 'Color', colors(iphs,:).^(1/fi)); hold on; end

    nexttile(3);
    semilogy(f(fplti,:), squeeze(dsc(1,2,:)), 'Color', colors(1,:).^(1/fi)); hold on; 
    semilogy(f(fplti,:), squeeze(dsc(1,3,:)), 'Color', colors(2,:).^(1/fi)); hold on; 
    semilogy(f(fplti,:), squeeze(dsc(2,3,:)), 'Color', colors(3,:).^(1/fi)); hold on; 


end


% finish up figure formatting
figure(figXf);
for iphs = 1:NPHS
    nexttile(iphs); ylabel('Connectivity'); title(['phase ' num2str(iphs)]);
end
legend(num2str((1:NPHS)'), Location='east'); 
xlabel(['phase ' num2str(fplti)]);


figure(figKC);
nexttile(1); title('$K_v/\phi$'); 
nexttile(2); title('$K_f/\phi$');
nexttile(3); title('$\phi^2/C_v$'); xlabel(['phase ' num2str(fplti)]);
nexttile(4); title('$\phi^2/C_f$'); xlabel(['phase ' num2str(fplti)]);
legend(num2str((1:NPHS)'), Location='east'); 


figure(figWgt);
nexttile(1); title('Pressure weights'); legend(num2str((1:NPHS)'), Location='east'); 
nexttile(2); title('Velocity weights');
nexttile(3); title('seg-comp length') ; xlabel(['phase ' num2str(fplti)]);
legend('2-1', '3-1', '3-2');


end