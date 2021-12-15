% make plots for permission functions and segregation-compaction lengths
% for olv-bas system, testing different d0 for basalt
% YQW, 13 May 2021

clear all;
addpath(genpath('../utils/'));

%% material properties

NPHS = 2;                   % number of phases

PHS  = {'olv','bas'};
rho0 = [ 3200; 2700];       % pure-phase densities
eta0 = [1e+18;1e+02];       % pure-phase viscosities
d0   = [ 5e-3; 5e-3];       % characteristic size of local-scale phase constituents

% set permission weight parameters for coefficient closure model
A = [ 0.6945, 0.1832; 0.5360, 0.1834;];  % permission slopes
B = [ 0.6906, 0.3094; 0.9993, 0.0007;];  % permission step locations
C = [ 0.6889, 0.1750; 0.8154, 1.5642;];  % permission step widths

%% initialize phase fractions

N    = 1e5;
fmid = round(0.5*N);

folv = linspace(0,1,N);
fbas = 1 - folv;
f    = [folv; fbas];

%% calculate permissions for different d0 basalt

d0bas = 1e-3*[0.1,0.5,1,5,10];
Nd0   = length(d0bas);

omvc  = zeros(NPHS, N, Nd0); omfc = omvc;
dsc   = zeros(NPHS, NPHS, N, Nd0);


for di = 1:Nd0
    d0(2) = d0bas(di);
    
    [dsc(:,:,:,di), Kv, Kf, Cv, Cf, Xf] = SegCompLength(f, eta0, d0, A, B, C);
    omvc(:,:,di) = Cv./sum(Cv);
    omfc(:,:,di) = Cf./sum(Cf);
    
end

%% plot

figure;
colors = parula(Nd0+1);
set(gcf,'defaultaxescolororder',repmat(colors(1:end-1,:),2,1));
set(gcf,'Position',[400,200,500,800]);
hAx = tight_subplot(3,1,0.05,[0.08,0.02],[0.15,0.05]);

axes(hAx(1));
semilogy(fbas, squeeze(dsc(1,2,:,:)));
hold on; plot(xlim, d0(1)*ones(1,2), 'k:'); hold off;
hleg = legend(num2str(1e3*d0bas(:)), 'location', 'northeast');
title(hleg, '$d_0^{bas}$ (mm)');
ylabel('Seg-comp length [m]');
ylim([1e-4,1e4]);
set(gca,'YTick',10.^(-4:2:4));


axes(hAx(2));
hsol = plot(fbas, squeeze(omvc(1,:,:))); hold on;
hliq = plot(fbas, squeeze(omvc(2,:,:)), '--'); hold off;
ylabel('Velocity weights');
legend([hsol(1), hliq(1)], {'solid','liquid'},'Location','west');

axes(hAx(3));
plot(fbas, squeeze(omfc(1,:,:))); hold on;
plot(fbas, squeeze(omfc(2,:,:)), '--'); hold off;
ylabel('Pressure weights');
xlabel('Liquid fraction $\phi^\ell$');

%% propagate effect of d0 into density calculations

rhoCf = sum(omfc.*rho0,1);
rhoCv = sum(omvc.*rho0,1);
rhobulk = sum(f.*rho0,1);

drho0 = rho0 - rhoCf;



%% plot densities

figure;
diplt = [1,3,4];
set(gcf,'Position',[500,500,1000,1000]);
t = tiledlayout(length(diplt),2,'TileSpacing','compact','Padding','compact');
colors = lines(3);
lbl1 = {'HorizontalAlignment', 'left' , 'FontSize', 18};
lbl2 = {'HorizontalAlignment', 'right', 'FontSize', 20, 'VerticalAlignment', 'top','Rotation',15};

for di = diplt
    nexttile;
    h = plot(fbas, rhobulk, fbas, rhoCf(:,:,di), fbas, rhoCv(:,:,di));
    hold on;
    plot(xlim, rho0(1)*ones(1,2), 'k-'); 
    plot(xlim, rho0(2)*ones(1,2), 'k-'); 
    hold off;
    uistack(h,'top');
    ylim([min(rho0)-50, max(rho0)+50]);
    
    legend(h, {'$\rho_{bulk}$', '$\rho^*_{Cf}$', '$\rho^*_{Cv}$'}, 'Location', 'west');
    
    
    % line labels
    text(0.8,rho0(1),'olivine','VerticalAlignment','bottom',lbl1{:});
    text(0.8,rho0(2),'basalt','VerticalAlignment','top',lbl1{:});
    
    xlabel('Liquid fraction, $\phi^\ell$');
    ylabel('Density [kg/m$^3$]');
    title(['$d_0^{bas}=$' num2str(1e3*d0bas(di)), ' mm']);
    
    
    
    nexttile;
    plot(fbas, rho0(1) - rhobulk, 'Color', colors(1,:)); hold on;
    plot(fbas, rho0(2) - rhobulk, 'Color', colors(1,:));
    plot(fbas, rhoCf(:,:,di) - rhobulk, 'Color', colors(2,:));
    plot(fbas, rhoCv(:,:,di) - rhobulk, 'Color', colors(3,:)); 
    plot(xlim, [0,0], 'k-', 'LineWidth', 0.8);
    hold off;
    
    % line labels
    text(0.95, rho0(1) - rhobulk(floor(0.05*N)), '$\rho^{olv}  - \rho_{bulk}$','color',colors(1,:),lbl2{:});
    text(0.95, rho0(2) - rhobulk(floor(0.05*N)), '$\rho^{bas}  - \rho_{bulk}$','color',colors(1,:),lbl2{:});
    text(0.30, rho0(2) - rhobulk(floor(0.45*N)), '$\rho^*_{Cf} - \rho_{bulk}$','color',colors(2,:),lbl2{:});
    text(0.30, rho0(1) - rhobulk(floor(0.75*N)), '$\rho^*_{Cv} - \rho_{bulk}$','color',colors(3,:),lbl2{:});
    
    xlabel('Liquid fraction, $\phi^\ell$');
    ylabel('$\Delta$ density [kg/m$^3$]');
    title(['$d_0^{bas}=$' num2str(1e3*d0bas(di)), ' mm']);


end


title(t,'Density (left) and density differences (right)', 'fontsize', 20, 'interpreter', 'latex');












