% make plots for permission functions and segregation-compaction lengths
% for plg-dac-mvp system
% YQW, 19 March 2021

clear all; close all;
Addpaths;

%%  set pure phase properties

PHS  = {'plg','dac','mvp'}; % phase names
NPHS = length(PHS);

rho0 = [3000 ;2500; 200];   % pure-phase densities
eta0 = [1e+16;1e+2;1e-3];   % pure-phase viscosities
d0   = [5e-3 ;5e-3;5e-3];   % characteristic size of local-scale phase constituents

% set permission weight parameters for coefficient closure model
% original
A = [ 0.60, 0.25, 0.30; 0.20, 0.20, 0.20; 0.20, 0.20, 0.20; ];  % permission slopes
B = [ 0.30, 0.15, 0.55; 0.48, 0.02, 0.50; 0.80, 0.08, 0.12; ];  % permission step locations
C = [ 0.20, 0.20, 0.20; 0.60, 0.60, 0.12; 0.20, 0.25, 0.50; ];  % permission step widths

% inspired by olv-bas
% A = [ 0.69, 0.18, 0.30; 0.54, 0.18, 0.20; 0.20, 0.20, 0.20; ];  % permission slopes
% B = [ 0.55, 0.30, 0.15; 0.48, 0.02, 0.50; 0.80, 0.08, 0.12; ];  % permission step locations
% C = [ 0.10, 0.18, 0.20; 0.82, 0.40, 0.12; 0.20, 0.25, 0.50; ];  % permission step widths

% inspired by parmigiani 2017
A = [ 0.60, 0.25, 0.30; 0.20, 0.20, 0.20; 0.20, 0.20, 0.20; ];  % permission slopes
B = [ 0.30, 0.15, 0.55; 0.48, 0.02, 0.50; 0.52, 0.43, 0.05; ];  % permission step locations
C = [ 0.20, 0.20, 0.20; 0.60, 0.60, 0.12; 1.00, 0.25, 1.50; ];  % permission step widths


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

ax = Plot3PhasePerm(f, Xf, PHS);
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



%% plot dsc 2D

PlotPerm2D(3, [0.1,0.05,0.01], eta0, d0, A, B, C, 2);

