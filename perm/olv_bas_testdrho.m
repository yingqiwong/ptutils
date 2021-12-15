% make plots for permission functions and segregation-compaction lengths
% for olv-bas system
% YQW, 18 March 2021

clear all;

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

%% calculate permissions

% get diffusivity contrasts
kv = eta0;          % momentum diffusivity
kf = d0.^2./eta0;   % volume diffusivity
Mv = kv.'./kv;      % momentum diffusivity ratios
Mf = kf.'./kf;      % volume diffusivity ratios

% get permission weights
F  = permute(repmat(f,1,1,1,NPHS),[4,1,2,3]);
Sf = (F./B).^(1./C);  Sf = Sf./sum(Sf,2);
Xf = sum(A.*Sf,2).*F + (1-sum(A.*Sf,2)).*Sf;

% get momentum and volume permissions
thtv = squeeze(prod(Mv.^Xf,2));  
thtf = squeeze(prod(Mf.^Xf,2));  

% get momentum and volume flux and transfer coefficients
Kv =    f .*eta0       .*thtv;
Kf =    f .*d0.^2./eta0.*thtf;
Cv = (1-f)./d0.^2.*Kv;
Cf = (1-f)./d0.^2.*Kf;

%% calculate segregation-compaction length

% get segregtation-compaction length scales
delta0(1,:) = f(2,:).*f(1,:)./sqrt(Cv(2,:).*Cf(1,:));
delta0(2,:) = f(1,:).*f(2,:)./sqrt(Cv(1,:).*Cf(2,:));

%% calculate velocity and pressure scales

omfc  =  Cf./sum(Cf,1);
rhoCf = sum(omfc.*rho0,1);

omvc  =  Cv./sum(Cv,1);
rhoCv = sum(omvc.*rho0,1);

rhobulk = sum(f.*rho0,1);

drho0 = rho0 - rhoCf;
% drho0 = rho0 - rhobulk;
% drho0   = rhobulk - rhoCf;
% drho0 = min(abs(drho0), [], 1);

% drho0 = 500;
g     = 9.81;
D     = 1;

% velocity scales
% uRef  = -500.*sign(folv-0.5).*g.*D.^2./sum(Kv,1);
uRef(:,:,1)  = -500.*sign(folv-0.5).*g.*D.^2./sum(Kv,1);
uRef(:,:,2)  = -(rhobulk-rhoCf).*sign(folv-0.5).*g.*D.^2./sum(Kv,1);
uRef(:,:,3)  = -(rhobulk-rhoCv).*sign(folv-0.5).*g.*D.^2./sum(Kv,1);


% uSegr = -f.^2.*drho0.*g./Cv;
uSegr(:,:,1) = -f.^2.*(rho0-rhobulk ).*g./Cv;
uSegr(:,:,2) = -f.^2.*(rhobulk-rhoCf).*g./Cv;
uSegr(:,:,3) = -f.^2.*(rhobulk-rhoCv).*g./Cv;

% pressure scales
% pRef  = 500.*g.*D;
pRef(:,:,1) = 500*g*D*ones(1,N);
pRef(:,:,2) = (rhobulk-rhoCf)*g*D;
pRef(:,:,3) = (rhobulk-rhoCv)*g*D;


% pComp = f.^2.*drho0.*g.*D./Cf./sum(Kv,1); 
% pComp(:,:,1) = f.^2.*(rho0-rhobulk ).*g.*D./Cf./sum(Kv,1); 
% pComp(:,:,2) = f.^2.*(rhobulk-rhoCf).*g.*D./Cf./sum(Kv,1); 
% pComp(:,:,3) = f.^2.*(rhobulk-rhoCv).*g.*D./Cf./sum(Kv,1); 

pComp(:,:,1) = (rho0-rhobulk ).*g.*delta0;
pComp(:,:,2) = (rhobulk-rhoCf).*g.*delta0;
pComp(:,:,3) = (rhobulk-rhoCv).*g.*delta0;


%% plot densities

figure;
set(gcf,'Position',[500,500,1000,400]);

subplot(121);
plot(fbas, rhobulk, fbas, rhoCf, fbas, rhoCv);
hold on;
h = plot(xlim, rho0(1)*ones(1,2), 'k-'); uistack(h, 'bottom');
h = plot(xlim, rho0(2)*ones(1,2), 'k-'); uistack(h, 'bottom');
hold off;
ylim([min(rho0)-50, max(rho0)+50]);

% line labels
lbl = {'HorizontalAlignment', 'left', 'FontSize', 20};
text(0.8,rho0(1),'olivine','VerticalAlignment','bottom',lbl{:});
text(0.8,rho0(2),'basalt','VerticalAlignment','top',lbl{:});
text(0.52,rhoCv(fmid),'$\rho^*_{Cv}$',lbl{:});
text(0.52,rhoCf(fmid),'$\rho^*_{Cf}$',lbl{:});
text(0.18,rhobulk(floor(0.7*N)),'$\rho_{bulk}$',lbl{:});

xlabel('Liquid fraction, $\phi^\ell$');
ylabel('Density [kg/m$^3$]');
title('(a) Densities in system');

subplot(122);
colors = lines(3);
plot(fbas, rho0(1) - rhobulk, '-', 'Color', colors(1,:)); hold on;
plot(fbas, rho0(2) - rhobulk, '--', 'Color', colors(1,:)); 
plot(fbas, rhoCf - rhobulk, '-', 'Color', colors(2,:));
plot(fbas, rhoCv - rhobulk, '-', 'Color', colors(3,:)); hold off;

% line labels
text(0.55, rho0(1) - rhobulk(floor(0.34*N)), '$\rho^{olv} - \rho_{bulk}$','Rotation',25,lbl{:});
text(0.7, rho0(2) - rhobulk(floor(0.4*N)), '$\rho^{bas} - \rho_{bulk}$','Rotation',25,lbl{:});
text(0.2, rhoCf(floor(0.6*N)) - rhobulk(floor(0.6*N)), '$\rho^*_{Cf} - \rho_{bulk}$','Rotation',25,lbl{:});
text(0.25, 0, '$\rho^*_{Cv} - \rho_{bulk}$','Rotation',25,lbl{:});

xlabel('Liquid fraction, $\phi^\ell$');
ylabel('$\Delta$ density [kg/m$^3$]');
title('(b) Density differences');


%% plot coefficients

HdrStyle = {'Units','normalized','FontSize',20,'HorizontalAlignment','center'};
LabStyle = {'Units','normalized','FontSize',18,'VerticalAlignment' ,'bottom'};
TtlStyle = {'Units','normalized','FontSize',18};
LegStyle = {'Location','east'   ,'FontSize',16,'Box','off'};

figure;
set(gcf,'Position',[400,400,800,800]);
hAx = tight_subplot(3,2,[0.08,0.03],[0.06,0.04],0.1);

% plot permission weights
axes(hAx(1));
plot(fbas, squeeze(Xf(1,1,:)), fbas, squeeze(Xf(2,2,:)));
ht=title('(a) Intra-phase weights',TtlStyle{:});
legend('solid', 'liquid', LegStyle{:});
ylabel('weights, $X_\phi^{ii}$');
set(gca,'YTick', 0:0.2:1);

axes(hAx(2));
plot(fbas, squeeze(Xf(1,2,:)), fbas, squeeze(Xf(2,1,:)));
ylabel('weights, $X_\phi^{ik}$');
ht=title('(b) Inter-phase weights',TtlStyle{:});
set(gca,'YAxisLocation','right');
set(gca,'YTick', 0:0.2:1);

axes(hAx(3));
semilogy(fbas, Kv./f);
title('(c) Effective viscosity', TtlStyle{:});
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

%% plot segregation-compaction lengths and velocity scales

figure;
set(gcf,'Position', [500,500,900,350]);
hAx = tight_subplot(1,2,0.1,[0.15,0.08],0.08);

axes(hAx(1));
semilogy(fbas, delta0, 'k-');
hold on; plot(xlim, d0(2)*ones(1,2), 'k:'); hold off;
text(0.1, d0(2), 'grain size', 'FontSize', 16, 'VerticalAlignment', 'bottom');
ylim([1e-3, 1e4]);
xlabel('Liquid fraction, $\phi^\ell$');
ylabel('seg-comp length, $\delta_{sc}^{bas-olv}$ [m]');
ht = title('(a)', 'HorizontalAlignment', 'left'); ht.Position(1) = 0;


axes(hAx(2));
hseg(1) = semilogy(fbas, -uSegr(1,:), '--'); hold on;
hseg(2) = semilogy(fbas,  uSegr(2,:));

% Dvec = [1;100;1e4];
% Dlbl = strcat({'$\ell_0 = 10^'}, num2str(log10(Dvec), '%.0f'), '$ m');

Dfac = [0.1;1;10;100];
Dvec = Dfac.*delta0(1,:);
Dlbl = num2str(1./Dfac);

uRefD = uRef.*Dvec.^2;
% href = semilogy(fbas(1:fmid), uRefD(:,1:fmid), 'color', 0.8*ones(1,3));
% semilogy(fbas(fmid+1:end), -uRefD(:,fmid+1:end), '--', 'color', 0.8*ones(1,3));

pos = uRef>=0;
href = semilogy(fbas( pos),  uRefD(:, pos), '-' , 'color', 0.8*ones(1,3));
href = semilogy(fbas(~pos), -uRefD(:,~pos), '--', 'color', 0.8*ones(1,3));

fbi = find(fbas<0.25,1);
text(fbas(fbi)*ones(size(Dfac)), -uRefD(:,fbi), Dlbl, 'FontSize',18, 'VerticalAlignment', 'middle');
hold off;
uistack(hseg, 'top');
xlabel('Liquid fraction, $\phi^\ell$');
ylabel('velocity magnitude [m/s]');
% ylim([1e-20,1e0]);
legend([href(1), hseg], 'reference, $u^*$', 'solid segr, $u_{segr}^s$', 'liquid segr, $u_{segr}^\ell$', 'Location', 'southeast', 'box', 'off');
ht = title('(b)', 'HorizontalAlignment', 'left'); ht.Position(1) = 0;

SaveFigure('Figures/olvbas_dscvelscale_rhoCv');


%% plot velocity and pressure scales

figure;
set(gcf,'Position',[500,500,900,600]);

hAx = tight_subplot(2,2,0.1,0.1,0.1);
lspec = {'LineStyle', 'none', 'MarkerSize', 12, 'linewidth', 3.5};

colororder(repmat(lines(3),2,1));

axes(hAx(1));
hold on;
plot(fbas, squeeze(abs(uRef).*100.*delta0(1,:).^2));
hold off;
% ylim([1e-8,1e-2]);
set(gca,'YScale','log','YTickLabelMode','auto','XTickLabelMode','auto');
ylabel('max $ |w^*|$ [m/s]');
title('(a) Vertical reference velocity');
hleg = legend('500','$\rho^*_{Cf} - \rho_{bulk}$','$\rho^*_{Cv} - \rho_{bulk}$','location','southeast','box','off');
title(hleg, '$\Delta \rho_0$');


axes(hAx(2));
plot(fbas, squeeze(abs(uSegr(1,:,:))), '--');
hold on;
plot(fbas, squeeze(abs(uSegr(2,:,:))), '-');
hold off;
% ylim([1e-20,1e0]); set(gca,'YTick',10.^(-20:5:0));
set(gca,'YScale','log','YTickLabelMode','auto','XTickLabelMode','auto');
hold off;
hleg = legend('$\rho^i - \rho_{bulk}$','$\rho^*_{Cf} - \rho_{bulk}$','$\rho^*_{Cv} - \rho_{bulk}$','location','southeast','box','off');
title(hleg, '$\Delta \rho_0$');
ylabel('max $ |w^i_\Delta|$ [m/s]');
title('(b) Vertical segregation velocity');


axes(hAx(3));
hold on;
semilogy(fbas, squeeze(abs(pRef*100.*delta0(1,:))));
hold off;
ylim([1e0,1e10]);
set(gca,'YScale','log','YTickLabelMode','auto','XTickLabelMode','auto');
set(gca,'YTick',10.^(-2:2:8));
xlabel('Liquid fraction, $\phi^\ell$');
ylabel('max $ |p^*|$ [Pa]');
title('(c) Reference pressure');


axes(hAx(4));
hold on;
semilogy(fbas, squeeze(abs(pComp(1,:,:))) , '--');
semilogy(fbas, squeeze(abs(pComp(2,:,:))) , '-');
hold off;
% ylim([1e-8,1e6]);
set(gca,'YScale','log','YTickLabelMode','auto','XTickLabelMode','auto');
set(gca,'YTick',10.^(-10:2:10));
xlabel('Liquid fraction, $\phi^\ell$');
ylabel('max $ |p_\Delta^i|$ [Pa]');
title('(d) Compaction pressure');


sgtitle('D = 100 $\delta_{sc}^{s\ell}$','FontSize',20);





























