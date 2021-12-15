function [uRef, uSegr, pRef, pComp, drho0] = CalcVelPressureScales (f, D, Kv, Cv, Cf, rho0)
% 
% [uRef, uSegr, pRef, pComp, pCompufac, drho0] = CalcVelPressureScales (f, D, Kv, Cv, Cf, rho0)
% 
% Calculates the expected scales of reference and segregation velocity, and
% reference and compaction pressures, based on scaling analysis from
% equations. 
% Still unresolved: what is the reference density?
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
% uRef      Reference velocity scale [NPHS x N or 1 x N]
% uSegr     Segregation velocity scale [NPHS x N]
% pRef      Reference pressure scale [NPHS x N or 1 x N]
% pComp     Compaction pressure scale [NPHS x N]
% pCompufac prefactor for compaction pressure scale [NPHS x N]
% drho0     Density scale used [NPHS x 1 or NPHS x N]
% 
% YQW, 19 April 2021

% calculate velocity-weighted density
omvc  =  Cv./sum(Cv,1);
rhoCv = sum(omvc.*rho0,1);

% pressure-weighted density
omfc  =  Cf./sum(Cf,1);
rhoCf = sum(omfc.*rho0,1);

rhobulk = sum(f.*rho0,1);

% Still unresolved: what is the reference density?
% drho0 = rhobulk - rhoCv;
% drho0 = rhobulk - rhoCf;
drho0 = (rho0 - rhobulk);
% drho0 = rho0 - rhoCv;
% drho0 = rho0(2) - rhobulk;
% drho0 = rho0 - rhobulk; 
% drho0 = max(rho0-rho0', [], 'all');
% drho0 = 500;

g     = 9.81;

% velocity scales
uRef  = max(abs(drho0)).*g.*D.^2./sum(Kv,1);
uSegr = -f.^2.*drho0.*g./Cv;

% pressure scales
pRef  = max(abs(drho0)).*g.*D;
pComp = f.^2.*drho0.*g.*D./Cf./sum(Kv,1); 
% pComp = f.^2.*(rho0-rhobulk).*g.*D./Cf./sum(Kv,1); 
end

