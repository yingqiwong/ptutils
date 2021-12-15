% make plots for permission functions and segregation-compaction lengths
% for plg-dac-mvp system
% YQW, 19 March 2021

addpath('../utils/');
clear variables;

%%  set pure phase properties (original)

PHS  = {'olv','plg','bas'}; % phase names
NPHS = length(PHS); 

rho0 = [ 3200; 2400; 2700]; % pure-phase densities
eta0 = [1e+18;1e+15;1e+02]; % pure-phase viscosities
d0   = [ 5e-3; 5e-3; 5e-3]; % characteristic size of local-scale phase constituents

% set permission weight parameters for coefficient closure model
A  =  [ 0.60, 0.30, 0.25; ...
        0.40, 0.25, 0.20; ...
        0.30, 0.25, 0.20; ];  % permission slopes
    
B  =  [ 0.45, 0.35, 0.20; ...
        0.35, 0.30, 0.35; ...
        0.45, 0.54, 0.01; ];  % permission step locations
    
C  =  [ 0.40, 0.40, 0.20; ...
        0.40, 0.40, 0.20; ...
        0.60, 0.20, 0.60; ];  % permission step widths
    
%% changed order of phases

% PHS  = {'olv','bas','plg'}; % phase names
% NPHS = length(PHS); 
% 
% rho0 = [ 3200; 2700; 2400]; % pure-phase densities
% eta0 = [1e+18;1e+02;1e+15]; % pure-phase viscosities
% d0   = [ 5e-3; 5e-3; 5e-3]; % characteristic size of local-scale phase constituents
% 
% % set permission weight parameters for coefficient closure model
% A  =  [ 0.60, 0.25, 0.30; ...
%         0.30, 0.20, 0.25; ...
%         0.40, 0.20, 0.25; ];  % permission slopes
%     
% B  =  [ 0.45, 0.20, 0.35; ...
%         0.45, 0.01, 0.54; ...
%         0.35, 0.35, 0.30; ];  % permission step locations
%     
% C  =  [ 0.40, 0.20, 0.40; ...
%         0.60, 0.60, 0.20; ...
%         0.40, 0.20, 0.40; ];  % permission step widths

%% alternative calibrations (for testing)

% PHS  = {'olv','plg','bas'}; % phase names
% NPHS = length(PHS); 
% 
% rho0 = [ 3200; 2400; 2700]; % pure-phase densities
% eta0 = [1e+18;1e+15;1e+02]; % pure-phase viscosities
% d0   = [ 5e-3; 5e-3; 5e-3]; % characteristic size of local-scale phase constituents
% 
% set permission weight parameters for coefficient closure model
% A  =  [ 0.60, 0.30, 0.25; ...
%         0.40, 0.25, 0.20; ...
%         0.30, 0.25, 0.20; ];  % permission slopes
%     
% B  =  [ 0.20, 0.05, 0.75; ...
%         0.35, 0.30, 0.35; ...
%         0.45, 0.54, 0.01; ];  % permission step locations
%     
% C  =  [ 0.40, 0.40, 0.20; ...
%         0.40, 0.40, 0.20; ...
%         0.60, 0.20, 0.60; ];  % permission step widths

%% initialize phase fractions

np    =  100;
f1  =  linspace(0,1,np);
f3  =  linspace(0,1,np);
[f3,f1]  =  meshgrid(f3,f1);
f3  =  f3(:);  f1  =  f1(:);  f2  =  1-f3-f1;

f               =  [f1,f2,f3].';
f(:,f(2,:)<0) = nan;

%% velocity and pressure scales

[dsc, Kv, Kf, Cv, Cf, Xf] = SegCompLength(f, eta0, d0, A, B, C);
[uRef, uSegr, pRef, pComp, pCompufac] = CalcVelPressureScales(f, 1, Kv, Cv, Cf, rho0);

%% plot connectivities

Plot3PhasePerm(f, Xf, PHS);
% SaveFigure('Figures/olvplgbas_varorderflip_calib');

%% plot permissions

kv = eta0;          % momentum diffusivity
Mv = kv.'./kv;      % momentum diffusivity ratios
thtv = squeeze(prod(Mv.^Xf,2));
Plot3PhaseCoeff(f, thtv, 'scl', 'log', 'PHS', PHS, 'cfname', {'\theta_v'});

Plot3PhaseCoeff(f, cat(3,Kv,Kf), 'scl', 'log', 'PHS', PHS, 'cfname', {'K_v','K_f'});
% SaveFigure('Figures/olvplgbas_varorderflip_K');

Plot3PhaseCoeff(f, cat(3,Cv,Cf), 'scl', 'log', 'PHS', PHS, 'cfname', {'C_v','C_f'});
% SaveFigure('Figures/olvplgbas_varorderflip_C');

%% plot weights

omCv = Cv./sum(Cv,1);
omCf = Cf./sum(Cf,1);
Plot3PhaseCoeff(f, cat(3,omCv,omCf), 'PHS', PHS, 'cfname', {'\omega_{Cv}','\omega_{Cf}'}, 'cflim', [0;1].*ones(1,2));
% SaveFigure('Figures/olvplgbas_varorderflip_omega');

%% plot seg-comp length

dscmat = [squeeze(dsc(1,2,:))'; squeeze(dsc(1,3,:))'; squeeze(dsc(2,3,:))'];
% plot3phasecoeff(f, dscmat, 'scl', 'log', 'PHS', PHS, 'cfname', {'\delta_{sc}'});

%% some 2d plots

%{
figure; 
set(gcf,'Position',[500,500,700,500]);
hAx = tight_subplot(2,2);

f3unique = unique(f3);
plgval = f3unique(100);

f1unique = unique(f1);
olvval = f1unique(1);
 
axes(hAx(1));
plot(f(2,f(3,:)==plgval), omCv(:,f(3,:)==plgval))
ylabel('Velocity weights');
title(['plg = ' num2str(plgval)]);
legend(PHS{:}, 'location', 'best', 'box', 'off');

axes(hAx(2));
plot(f(2,f(1,:)==olvval), omCv(:,f(1,:)==olvval))
title(['olv = ' num2str(olvval)]);

axes(hAx(3));
plot(f(2,f(3,:)==plgval), omCf(:,f(3,:)==plgval))
ylabel('Pressure weights'); xlabel('basalt fraction');

axes(hAx(4));
plot(f(2,f(1,:)==olvval), omCf(:,f(1,:)==olvval))
xlabel('basalt fraction');
%}






