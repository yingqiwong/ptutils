function [dzi, ft] = CrossCorrPhaseFracs (t, x, f, ti, xi)
%
% [t, x, ft, dxi] = CrossCorrPhaseFracs (t, x, f, ti)
% cross-correlate vertical phase fraction profiles to see how far the 
% profile has moved over a few timesteps
% 
% example: CrossCorrPhaseFracs(t, x, f, 20:30);
%
% INPUTS
% t         simulation times [Nf x 1]
% x         x positions [1 x N]
% f         phase fractions [NPHS x N x N x Nf]
% ti        time indices to calculate cross correlations [Nti x 1]
%
% OUTPUTS
% dzi       grid cell offset at each time step [Nti x 1]
% ft        translated phase fraction profiles [NPHS x N x N x Nti]
%
% DEFAULT OPTIONS
%
% YQW, 17 June 2021


% get sizes
Nt = length(ti);
N  = size(f,2);

% get profiles
t = t(ti);
f = f(:,:,:,ti);

% perform cross correlations along center profile to get offsets
if nargin<5, xi = floor(0.5*N); end
dzi = zeros(size(t));

for j = 2:length(t)
    r = crosscorrelation(f(1,:,xi,j), f(1,:,xi,1));
    [~,dzi(j)] = max(r); 
end


ft = f;
% now translate profiles
for j = 2:length(t)
    zind  = [(dzi(j)):N, 1:(dzi(j)-1)];
    zind(zind<=0) = [];
    ft(:,:,:,j) = f(:,zind,:,j);
end

% plot profiles
figure;
set(gcf,'defaultaxescolororder', copper(Nt));

subplot(121);
plot(squeeze(f(1,:,xi,:)), 1:N); hold on;
title('original profiles');
leg = legend(num2str(ti'), 'location', 'best');
title(leg, '$t_i$');

subplot(122);
plot(squeeze(ft(1,:,xi,:)), 1:N); hold on;
title('translated profiles');




end

function [r] = crosscorrelation (f1, f2)
% calculate cross-correlation between two vectors

r = zeros(size(f1));
N = length(f1);

for i = 0:(N-1)
    r(i+1) = 1/N*sum(f1.*f2);
    f2 = f2([end (1:end-1)]);
end

end

