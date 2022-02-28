% make plots for permission functions and segregation-compaction lengths
% for plg-dac-mvp system
% YQW, 19 March 2021

clear variables; close all;
Addpaths;

%%  set pure phase properties  

PHS  = {'plg','dac','fep'}; % phase names
NPHS = length(PHS); 


rho0 = [3000 ;2500;4000];   % pure-phase densities
eta0 = [1e+16;1e+3;1e+0];   % pure-phase viscosities
d0   = [1e-2 ;1e-2;1e-2];   % characteristic size of local-scale phase constituents

% set permission weight parameters for coefficient closure model
A = [ 0.60, 0.25, 0.30; 0.20, 0.20, 0.20; 0.20, 0.20, 0.20; ];  % permission slopes
B = [ 0.30, 0.15, 0.55; 0.48, 0.02, 0.50; 0.80, 0.08, 0.12; ];  % permission step locations
C = [ 0.20, 0.20, 0.20; 0.60, 0.60, 0.12; 0.20, 0.25, 0.50; ];  % permission step widths

%% initialize phase fractions
np    =  200;
f1  =  linspace(0,1,np);
f2  =  linspace(0,1,np);
[f2,f1]  =  meshgrid(f2,f1);
f2  =  f2(:);  f1  =  f1(:);  f3  =  1-f2-f1;

f               =  [f1,f2,f3].';
f(:,f(3,:)<0) = nan;

%% velocity and pressure scales

[dsc, Kv, Kf, Cv, Cf, Xf] = SegCompLength(f, eta0, d0, A, B, C);
[uRef, uSegr, pRef, pComp] = CalcVelPressureScales(f, 1, Kv, Cv, Cf, rho0);

%% plot permissions

Plot3PhasePerm(f, Xf, PHS);
% SaveFigure('Figures/plgdacfep_calib');

Plot3PhaseCoeff(f, cat(3,Kv,Kf), 'scl', 'log', 'PHS', PHS, 'cfname', {'K_v','K_f'});
% SaveFigure('Figures/plgdacfep_K');

Plot3PhaseCoeff(f, cat(3,Cv,Cf), 'scl', 'log', 'PHS', PHS, 'cfname', {'C_v','C_f'});
% SaveFigure('Figures/plgdacfep_C');

omCv = Cv./sum(Cv,1);
omCf = Cf./sum(Cf,1);
plot3phasecoeff(f, cat(3,omCv,omCf), 'PHS', PHS, 'cfname', {'\omega_{Cv}','\omega_{Cf}'}, 'cflim', [0;1].*ones(1,2));
% SaveFigure('Figures/plgdacfep_omega');

dscmat = [squeeze(dsc(1,2,:))'; squeeze(dsc(1,3,:))'; squeeze(dsc(2,3,:))'];
plot3phasecoeff(f, dscmat, 'scl', 'log', 'PHS', PHS, 'cfname', {'\delta_{sc}'});
