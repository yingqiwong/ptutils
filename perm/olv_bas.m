% make plots for permission functions and segregation-compaction lengths
% for olv-bas system
% YQW, 18 March 2021

clear all;
Addpaths

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

N    = 1e6;
fmid = round(0.5*N);

folv = linspace(0,1,N);
fbas = 1 - folv;
f    = [folv; fbas];

%% calculate permissions

[dsc, Kv, Kf, Cv, Cf, Xf] = SegCompLength(f, eta0, d0, A, B, C);
[uRef, uSegr, pRef, pComp, pCompufac] = CalcVelPressureScales(f, 1, Kv, Cv, Cf, rho0);
omvc = Cv./sum(Cv);
omfc = Cf./sum(Cf);

Kvbar = sum(Kv,1);
segnum = f.^2.*Kvbar./Cv;
compnum = f.^2./Cf./Kvbar;

%%
HdrStyle = {'Units','normalized','FontSize',20,'HorizontalAlignment','center'};
LabStyle = {'Units','normalized','FontSize',18,'VerticalAlignment' ,'bottom'};
TtlStyle = {'Units','normalized','FontSize',18};
LegStyle = {'Location','east'   ,'FontSize',16,'Box','off'};


%% just plot seg comp length

figure;
set(gcf,'Position', [500,500,500,350]);
semilogy(fbas, squeeze(dsc(1,2,:)), 'k-');
hold on; 
plot(xlim, d0(1)*ones(1,2), 'k:'); 
semilogy(fbas, 100*squeeze(dsc(1,2,:)), 'k--');
hold off;
text(0.1, d0(1), 'grain size', 'FontSize', 16, 'VerticalAlignment', 'bottom');
ylim([1e-3, 1e6]);
set(gca,'XTick', 0:0.2:1);
set(gca,'YTick', 10.^(-4:2:6),'YMinorTick','off');
xlabel('Liquid fraction, $\phi^\ell$');
ylabel('$\delta_{sc}^{bas-olv}$ [m]');
title('basalt seg - olivine comp length', TtlStyle{:});


%% plot coefficients

figure;
set(gcf,'Position',[400,400,800,800]);
hAx = tight_subplot(3,2,[0.08,0.03],[0.06,0.04],0.1);

% plot permission weights
axes(hAx(1));
plot(fbas, squeeze(Xf(1,1,:)), fbas, squeeze(Xf(2,2,:)));
ht=title('(a) Intra-phase weights',TtlStyle{:});
legend('solid-solid', 'liquid-liquid', LegStyle{:});
ylabel('weights, $X_\phi^{ii}$');
set(gca,'YTick', 0:0.2:1);

axes(hAx(2));
plot(fbas, squeeze(Xf(1,2,:)), fbas, squeeze(Xf(2,1,:)));
ylabel('weights, $X_\phi^{ik}$');
ht=title('(b) Inter-phase weights',TtlStyle{:});
legend('solid-liquid', 'liquid-solid', LegStyle{:});
set(gca,'YAxisLocation','right');
set(gca,'YTick', 0:0.2:1);

axes(hAx(3));
semilogy(fbas, Kv./f);
title('(c) Effective viscosity', TtlStyle{:});
legend('solid', 'liquid', LegStyle{:});
ylabel('$K_v^i/\phi^i$ [Pa s]');
ylim([1e2,1e18]);
set(gca,'YTick', 10.^(2:4:22),'YMinorTick','off');

% plot flux coeffs
axes(hAx(4));
semilogy(fbas, Kf./f);
title('(d) Volume diffusivity',TtlStyle{:}); 
ylabel('$K_\phi^i/\phi^i$ [m$^2$/Pa s]');
ylim([1e-24,1e-6]);
set(gca,'YTick',10.^(-24:4:-6),'YAxisLocation','right','YMinorTick','off');

axes(hAx(5));
semilogy(fbas, f.^2./Cv);
xlabel('Liquid fraction, $\phi^\ell$');
title('(e) Segregation coefficient', TtlStyle{:});
ylabel('$\phi^{i^2}/C_v^i$ [m$^2$/Pa s]');
ylim([1e-22,1e-2]);
set(gca,'YTick', 10.^(-22:4:-2),'YMinorTick','off');

% plot transfer coeffs
axes(hAx(6));
semilogy(fbas, f.^2./Cf);
xlabel('Liquid fraction, $\phi^\ell$');
title('(f) Compaction coefficient', TtlStyle{:});
ylabel('$\phi^{i^2}/C_\phi^i$ [Pa s]');
ylim([1e-2,1e22]);
set(gca,'YAxisLocation','right','YTick', 10.^(-2:4:22),'YMinorTick','off');

% ylim([1e-21,1e-3]);
% set(gca,'YTick',[1e-20,1e-15,1e-10,1e-5]);

% SaveFigure('Figures/olvbas_coeffs');

%% plot coefficients different way

figure;
set(gcf,'Position',[400,400,1200,550]);
hAx = tight_subplot(2,4,[0.12,0.06],[0.1,0.08],[0.05,0.03]);

% plot permission weights
axes(hAx(1));
plot(fbas, squeeze(Xf(1,1,:)), fbas, squeeze(Xf(2,2,:)));
ht=title('(a) Intra-phase weights',TtlStyle{:});
legend('solid-solid', 'liquid-liquid', LegStyle{:});
ylabel('weights, $X_\phi^{ii}$');
set(gca, 'YTick', 0:0.2:1, 'XTick',  0:0.2:1);

axes(hAx(5));
plot(fbas, squeeze(Xf(1,2,:)), fbas, squeeze(Xf(2,1,:)));
ylabel('weights, $X_\phi^{ik}$');
ht=title('(b) Inter-phase weights',TtlStyle{:});
legend('solid-liquid', 'liquid-solid', LegStyle{:});
set(gca, 'YTick', 0:0.2:1, 'XTick',  0:0.2:1);
xlabel('Liquid fraction, $\phi^\ell$');

axes(hAx(2));
semilogy(fbas, Kv);
title('(c) Momentum flux coefficient', TtlStyle{:});
legend('solid', 'liquid', LegStyle{:});
ylabel('$K_v^i$ [Pa s]');
set(gca, 'XTick',  0:0.2:1);
% ylim([1e2,1e18]);
% set(gca,'YTick', 10.^(2:4:22),'YMinorTick','off');

% plot flux coeffs
axes(hAx(6));
semilogy(fbas, Kf);
title('(d) Volume flux coefficient',TtlStyle{:}); 
ylabel('$K_\phi^i$ [m$^2$/Pa s]');
xlabel('Liquid fraction, $\phi^\ell$');
set(gca, 'XTick',  0:0.2:1);
ylim([1e-25,1e-5]);
set(gca,'YTick',10.^(-25:5:5));

axes(hAx(3));
semilogy(fbas, Cv);
title('(e) Momentum transfer coefficient', TtlStyle{:});
ylabel('$C_v^i$ [Pa s/m$^2$]');
set(gca, 'XTick',  0:0.2:1);
% ylim([1e-22,1e-2]);
% set(gca,'YTick', 10.^(-22:4:-2),'YMinorTick','off');

% plot transfer coeffs
axes(hAx(7));
semilogy(fbas, Cf);
xlabel('Liquid fraction, $\phi^\ell$');
title('(f) Volume transfer coefficient', TtlStyle{:});
ylabel('$C_\phi^i$ [/Pa s]');
xlabel('Liquid fraction, $\phi^\ell$');
set(gca, 'XTick',  0:0.2:1);
% ylim([1e-2,1e22]);
% set(gca,'YAxisLocation','right','YTick', 10.^(-2:4:22),'YMinorTick','off');

axes(hAx(4));
plot(fbas, Cv./sum(Cv,1));
set(gca, 'XTick',  0:0.2:1);
title('(g) Velocity weights', TtlStyle{:});
ylabel('$\omega_{C_\phi}^i$');

axes(hAx(8));
plot(fbas, Cf./sum(Cf,1));
xlabel('Liquid fraction, $\phi^\ell$');
title('(h) Pressure weights', TtlStyle{:});
ylabel('$\omega_{C_v}^i$');
xlabel('Liquid fraction, $\phi^\ell$');
set(gca, 'XTick',  0:0.2:1);

% ylim([1e-21,1e-3]);
% set(gca,'YTick',[1e-20,1e-15,1e-10,1e-5]);

SaveFigure('Figures/olvbas_coeffs');

%% plot coefficients

HdrStyle = {'Units','normalized','FontSize',20,'HorizontalAlignment','center'};
LabStyle = {'Units','normalized','FontSize',18,'VerticalAlignment' ,'bottom'};
TtlStyle = {'Units','normalized','FontSize',18};
LegStyle = {'Location','east'   ,'FontSize',16,'Box','off'};

figure;
set(gcf,'Position',[400,400,800,300]);
hAx = tight_subplot(1,2,[0.08,0.03],[0.17,0.1],0.1);

colors = lines(2);

% plot permission weights
axes(hAx(1));
plot(fbas, squeeze(Xf(1,1,:)), '-',  'Color', colors(1,:)); hold on;
plot(fbas, squeeze(Xf(1,2,:)), '--', 'Color', colors(1,:));
ht=title('Connectivity of the solid',TtlStyle{:});
legend('solid-solid', 'solid-liquid', LegStyle{:});
ylabel('weights, $X_\phi^{ik}$');
set(gca,'YTick', 0:0.2:1);
xlabel('Liquid fraction, $\phi^\ell$');

axes(hAx(2));
plot(fbas, squeeze(Xf(2,2,:)), '-',  'Color', colors(2,:)); hold on;
plot(fbas, squeeze(Xf(2,1,:)), '--', 'Color', colors(2,:));
ylabel('weights, $X_\phi^{ik}$');
ht=title('Connectivity of the liquid',TtlStyle{:});
legend('liquid-liquid', 'liquid-solid', LegStyle{:});
set(gca,'YAxisLocation','right');
set(gca,'YTick', 0:0.2:1);
xlabel('Liquid fraction, $\phi^\ell$');

% SaveFigure('Figures/olvbas_connectivity');

%% plot viscosity with an exponential model overlaid

fbas1 = 0.25;
fi = find(f(2,:)<fbas1,1);
lambda = [27;80;130];
eta1 = Kv(1,fi)./f(1,fi).*exp(lambda*fbas1);
eta = eta1.*exp(-lambda.*f(2,:));

figure;
set(gcf,'Position',[500,500,500,300]);
semilogy(fbas, Kv./f);
hold on; 
plot(0.3*ones(1,2), ylim, 'k-', 'linewidth', 0.5);
plot(fbas, eta, 'k:'); 
hold off;
text(0.8,2e6,num2str(lambda(1)),LabStyle{3:4});
text(0.5,2e6,num2str(lambda(2)),LabStyle{3:4});
text(0.33,2e6,num2str(lambda(3)),LabStyle{3:4});
title('Effective viscosity', TtlStyle{:});
legend('solid', 'liquid', LegStyle{:});
xlabel('Liquid fraction, $\phi^\ell$');
ylabel('$K_v^i/\phi^i$ [Pa s]');
ylim([1e2,1e18]);
set(gca,'YTick', 10.^(2:4:22),'YMinorTick','off');


%% plot segregation and compaction numbers

figure;
subplot(211); semilogy(fbas,  segnum); title('segregation number');
subplot(212); semilogy(fbas, compnum); title('compaction number');


%% plot segregation-compaction lengths and velocity scales

LegStyle = {'Location','west','FontSize',16,'Box','off'};

figure;
set(gcf,'Position', [500,500,1200,500]);
hAx = tight_subplot(2,3,0.08,[0.08,0.04],0.08);

hAx(1).Visible = 'off';

axes(hAx(2));
plot(fbas, omvc);
title('(a) Velocity weights', TtlStyle{:});
ylabel('$\omega_{C_v}$'); 
legend('solid', 'liquid', LegStyle{:});

axes(hAx(3));
plot(fbas, omfc);
title('(a) Pressure weights', TtlStyle{:});
ylabel('$\omega_{C_\phi}$'); 
legend('solid', 'liquid', LegStyle{:});


axes(hAx(4));
semilogy(fbas, squeeze(dsc(1,2,:)), 'k-');
hold on; plot(xlim, d0(1)*ones(1,2), 'k:'); hold off;
text(0.1, d0(1), 'grain size', 'FontSize', 16, 'VerticalAlignment', 'bottom');
ylim([1e-3, 1e4]);
set(gca,'YTick', 10.^(-4:2:4),'YMinorTick','off');
xlabel('Liquid fraction, $\phi^\ell$');
ylabel('$\delta_{sc}^{bas-olv}$ [m]');
title('(c) Segregation-compaction length', TtlStyle{:});


% assign domain sizes
Dfac   = [0.1;1;10;1000];
Dvec   = Dfac.*squeeze(dsc(1,2,:))';
Dlbl   = num2str(Dfac);
uRefD  = uRef.*Dvec.^2;
pRefD  = pRef.*Dvec;

axes(hAx(5));
hseg(1) = semilogy(fbas, abs(uSegr(1,:))); hold on;
hseg(2) = semilogy(fbas, abs(uSegr(2,:)));
href    = semilogy(fbas, abs(uRefD), '-' , 'color', 0.8*ones(1,3));

fbi = find(fbas<0.25,1);
text(fbas(fbi)*ones(size(Dfac)), -uRefD(:,fbi), Dlbl, 'FontSize',18, 'VerticalAlignment', 'middle');
hold off;
uistack(hseg, 'top');
xlabel('Liquid fraction, $\phi^\ell$');
ylabel('magnitude [m/s]');
% ylim([1e-20,1e0]);
legend([href(1), hseg], 'reference, $u^*$', 'solid segr, $u_{segr}^s$', 'liquid segr, $u_{segr}^\ell$', 'Location', 'southeast', 'box', 'off');
title('(d) Velocity scales', TtlStyle{:});


axes(hAx(6));
hseg(1) = semilogy(fbas, abs(pComp(1,:).*Dvec(3,:))); hold on;
hseg(2) = semilogy(fbas, abs(pComp(2,:).*Dvec(3,:)));
href    = semilogy(fbas, abs(pRefD), '-' , 'color', 0.8*ones(1,3));

fbi = find(fbas<0.25,1);
text(fbas(fbi)*ones(size(Dfac)), -pRefD(:,fbi), Dlbl, 'FontSize',18, 'VerticalAlignment', 'middle');
hold off;
uistack(hseg, 'top');
xlabel('Liquid fraction, $\phi^\ell$');
ylabel('magnitude [m/s]');
% ylim([1e-20,1e0]);
legend([href(1), hseg], 'reference, $u^*$', 'solid segr, $u_{segr}^s$', 'liquid segr, $u_{segr}^\ell$', 'Location', 'southeast', 'box', 'off');
title('(e) Pressure scales', TtlStyle{:});

% SaveFigure('Figures/olvbas_dscvelscale_rhoCv');



%% different way of plotting velocity scales

figure;
set(gcf,'Position', [500,500,500,700]);
hAx = tight_subplot(2,1,0.15,[0.08,0.05],[0.15,0.03]);

axes(hAx(1));
semilogy(fbas, squeeze(dsc(1,2,:)) , 'k-');
hold on; plot(xlim, d0(1)*ones(1,2), 'k:');
text(0.1, d0(1), 'solid grain size', 'FontSize', 16, 'VerticalAlignment', 'bottom');
ylim([1e-3, 1e4]);
set(gca,'YTick', 10.^(-4:2:4),'YMinorTick','off');
xlabel('Liquid fraction, $\phi^\ell$');
ylabel('$\delta_{sc}^{bas-olv}$ [m]');
title('basalt seg -- olivine comp length');

axes(hAx(2));
% assign domain sizes
Dfac   = [0.01;1;100];
Dvec   = Dfac.*squeeze(dsc(1,2,:))';
Dlbl   = num2str(Dfac);
uRefD  = uRef.*Dvec.^2;
pRefD  = pRef.*Dvec;
uRefFixedD  = uRef.*1e3;

% segregation velocities
hseg(1) = semilogy(fbas, abs(uSegr(1,:))); hold on;
hseg(2) = semilogy(fbas, abs(uSegr(2,:)));
txtstyle = {'FontSize',18, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center'};
text(0.8, 1e-2, '$u_\Delta^\ell$', txtstyle{:}, 'color', hseg(2).Color);
text(0.8, 4e-6, '$u_\Delta^s$'   , txtstyle{:}, 'color', hseg(1).Color);

% reference velocity by delta/D
% href = semilogy(fbas, abs(uRefD), '-' , 'color', 0.8*ones(1,3));
% fbi  = find(fbas<0.1,1);
% text(fbas(fbi)*ones(size(Dfac)), abs(uRefD(:,fbi)), Dlbl, txtstyle{:});
text(0.8, 100, 'ref $u^*$', txtstyle{:});

% reference velocity by absolute distance
semilogy(fbas, abs(uRefFixedD), 'k-');
fbi = find(fbas<0.4,1);
text(fbas(fbi), 3e3, 'D=1km' , txtstyle{:});


hold off;
uistack(hseg, 'top');
xlabel('Liquid fraction, $\phi^\ell$');
ylabel('magnitude [m/s]');
ylim([1e-20,1e5]);
title('velocity scales');

% SaveFigure('Figures/olvbas_dscvelscale_segreffixedD');


%% some sort of regime diagram?

Dvec = logspace(-1,6,101)';
Rsc = (squeeze(dsc(1,2,:))'./Dvec).^2;

figure; imagesc(fbas, log10(Dvec), log10(Rsc)); colorbar
set(gca,'ydir','normal');


