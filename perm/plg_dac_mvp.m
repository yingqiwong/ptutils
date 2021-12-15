% make plots for permission functions and segregation-compaction lengths
% for plg-dac-mvp system
% YQW, 19 March 2021

clear variables; close all;
addpath(genpath('../utils/'));

%%  set pure phase properties

PHS  = {'plg','dac','mvp'}; % phase names
NPHS = length(PHS);

rho0 = [3000 ;2500; 200];   % pure-phase densities
eta0 = [1e+16;1e+2;1e-3];   % pure-phase viscosities
d0   = [5e-3 ;5e-3;5e-3];   % characteristic size of local-scale phase constituents

% set permission weight parameters for coefficient closure model
% original
% A = [ 0.60, 0.25, 0.30; 0.20, 0.20, 0.20; 0.20, 0.20, 0.20; ];  % permission slopes
% B = [ 0.30, 0.15, 0.55; 0.48, 0.02, 0.50; 0.80, 0.08, 0.12; ];  % permission step locations
% C = [ 0.20, 0.20, 0.20; 0.60, 0.60, 0.12; 0.20, 0.25, 0.50; ];  % permission step widths

% inspired by olv-bas
A = [ 0.69, 0.18, 0.30; 0.54, 0.18, 0.20; 0.20, 0.20, 0.20; ];  % permission slopes
B = [ 0.55, 0.30, 0.15; 0.48, 0.02, 0.50; 0.80, 0.08, 0.12; ];  % permission step locations
C = [ 0.10, 0.18, 0.20; 0.82, 0.40, 0.12; 0.20, 0.25, 0.50; ];  % permission step widths

%% plot dsc 2D

f3 = 0.30;
N = 1001;

f2 = linspace(0, 1-f3, N);
f  = [1-f2-f3; f2; f3.*ones(1,N)];
f(f<0) = nan;
[dsc, Kv, Kf, Cv, Cf, Xf] = SegCompLength(f, eta0, d0, A, B, C);


figure;
set(gcf,'Position',[500,200,800,500]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
nexttile; semilogy(f2,    Kv./f); ylabel('$K_v/\phi$');
legend('solid','liquid','gas','Location','best');
nexttile; semilogy(f2,    Kf./f); ylabel('$K_f/\phi$');
nexttile; semilogy(f2, f.^2./Cv); ylabel('$\phi^2/C_v$');
nexttile; semilogy(f2, f.^2./Cf); ylabel('$\phi^2/C_f$');
sgtitle(['gas fraction = ' num2str(f3*100,'%.0f') '%']);
SaveFigure(['Figures/plgdacmvp_2D_coef_gas' num2str(f3*100,'%.0f') ]);


figure;
set(gcf,'Position',[500,200,500,900]);

subplot(311);
plot(f2, Cf./sum(Cf));
ylabel('pressure weights');
legend('solid','liquid','gas','Location','best');
title(['gas fraction = ' num2str(f3*100,'%.0f') '\%']);

subplot(312);
plot(f2, Cv./sum(Cv));
ylabel('velocity weights');

subplot(313);
semilogy(f2, squeeze(dsc(1,2,:)));      hold on;
semilogy(f2, squeeze(dsc(1,3,:)));
semilogy(f2, squeeze(dsc(2,3,:)));      hold off;
ylim([1e-8,1e4]);
ylabel('seg-comp length, $\delta_{sc}^{ik}$ [m]');
xlabel('Liquid fraction $\phi^\ell$');
leg = legend('liquid-solid','gas-solid','gas-liquid', 'box', 'off');

SaveFigure(['Figures/plgdacmvp_2D_dsc_gas' num2str(f3*100,'%.0f') ]);

%% initialize phase fractions

np    =  200;
f1  =  linspace(0,1,np);
f2  =  linspace(0,1,np);
[f2,f1]  =  meshgrid(f2,f1);
f2  =  f2(:);  f1  =  f1(:);  f3v  =  1-f2-f1;

f             =  [f1,f2,f3v].';
f(:,f(3,:)<0) = nan;

%% velocity and pressure scales

[dsc, Kv, Kf, Cv, Cf, Xf] = SegCompLength(f, eta0, d0, A, B, C);
[uRef, uSegr, pRef, pComp, pCompufac] = CalcVelPressureScales(f, 1, Kv, Cv, Cf, rho0);

%% plot permissions

Plot3PhasePerm(f, Xf, PHS);
% SaveFigure('Figures/plgdacmvp_calib');

Plot3PhaseCoeff(f, cat(3,Kv,Kf), 'scl', 'log', 'PHS', PHS, 'cfname', {'K_v','K_f'});
% SaveFigure('Figures/plgdacmvp_K');

Plot3PhaseCoeff(f, cat(3,Cv,Cf), 'scl', 'log', 'PHS', PHS, 'cfname', {'C_v','C_f'});
% SaveFigure('Figures/plgdacmvp_C');

omCv = Cv./sum(Cv,1);
omCf = Cf./sum(Cf,1);
Plot3PhaseCoeff(f, cat(3,omCv,omCf), 'PHS', PHS, 'cfname', {'\omega_{Cv}','\omega_{Cf}'}, 'cflim', [0;1].*ones(1,2));
% SaveFigure('Figures/plgdacmvp_omega');

dscmat = [squeeze(dsc(1,2,:))'; squeeze(dsc(1,3,:))'; squeeze(dsc(2,3,:))'];
Plot3PhaseCoeff(f, dscmat, 'scl', 'log', 'PHS', PHS, 'cfname', {'\delta_{sc}'});



