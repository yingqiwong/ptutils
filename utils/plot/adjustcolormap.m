function [cmapadj, v0, vlim] = adjustcolormap (cmap, v, v0)
%
% [cmap, v0, vlim] = adjustcolormap (cmap, v, v0, vlim)
%
% Adjust colormap colors so that the neutral color is at v0 even if
% vlim is not symmetric around v0.
% This function only makes sense if cmap is a diverging colorscale.
%
% INPUTS
% cmap      original colormap [Ncolors x 3]
% v         EITHER
%               a) vector of values represented by colormap [N x 1]
%               b) limits of values [2 x 1]
% v0        center of values [scalar]
%
% OUTPUTS
% cmapadj   adjusted colormap
% v0, vlim as in INPUTS
%
% YQW, 27 August 2021


% check definition of input v
if length(v) == 2
    vlim = v;
else
    vlim = [min(v), max(v)];
end


vdist = vlim - v0;
vdist = round(abs(vdist/vdist(1)));

cmid = floor(0.5*length(cmap));

if vdist(2)==vdist(1)
    % equal colors on each side
    cmapadj = cmap;
    
elseif vdist(2)>vdist(1)
    % more colors on upper side
    fac = vdist(2)/vdist(1);
    cmapadj = [cmap(1:fac:cmid,:); cmap(cmid+1:end,:)];
    
elseif vdist(2)<vdist(1)
    % more colors on lower side
    fac = vdist(1)/vdist(2);
    cmapadj = [cmap(1:cmid+1,:); cmap(cmid+2:fac:end,:)];
    
end




end