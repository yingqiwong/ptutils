function [uBar, uSegr, pBar, pComp, drhoseg] = CalcVelPressureScales (f, D, Kv, Cv, Cf, rho0)
% 
% [uRef, uSegr, pRef, pComp, pCompufac, drhoseg] = CalcVelPressureScales (f, D, Kv, Cv, Cf, rho0)
% 
% Calculates the expected scales of reference and segregation velocity, and
% reference and compaction pressures, based on scaling analysis from
% equations. 
% Be careful about the reference fields. They probably need to be rescaled
% after this using the phase fraction perturbations
% 
% INPUTS
% f         phase fractions [NPHS x N]
% D         domain size [1 x N]
% Kv        momentum flux coeff [NPHS x N]
% Cv        momentum transfer coeff [NPHS x N]
% Cf        volume transfer coeff [NPHS x N]
% rho0      phase density [NPHS x 1]
% 
% OUTPUTS
% uBar      Reference velocity scale [NPHS x N or 1 x N]
% uSegr     Segregation velocity scale [NPHS x N]
% pBar      Reference pressure scale [NPHS x N or 1 x N]
% pComp     Compaction pressure scale [NPHS x N]
% pCompufac prefactor for compaction pressure scale [NPHS x N]
% drhoseg   Density scale used for phase segregation [NPHS x 1 or NPHS x N]
% 
% YQW, 19 April 2021

% calculate velocity-weighted density
omvc  =  Cv./sum(Cv,1);
rhoCv = sum(omvc.*rho0,1);

% pressure-weighted density
omfc  =  Cf./sum(Cf,1);
rhoCf = sum(omfc.*rho0,1);

rhobulk = sum(f.*rho0,1);

drhobar = max(rho0-rho0',[],'all');
drhoseg = rho0 - rhobulk;

% gravity [m/s2]
g     = 9.81;

% velocity scales
uBar  = drhobar.*g.*D.^2./sum(Kv,1);
uSegr = -f.^2.*drhoseg.*g./Cv;

% pressure scales
pBar  = drhobar.*g.*D;
pComp = f.^2.*drhoseg.*g.*D./Cf./sum(Kv,1); 
end

