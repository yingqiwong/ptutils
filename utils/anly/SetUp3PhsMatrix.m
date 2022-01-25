function [f] = SetUp3PhsMatrix (np, ffixi, ffix)
% [f] = SetUp3PhsMatrix (np, ffixi, ffix)
% 
% 
% set up a matrix of phase fractions for THREE-PHASE mixtures. 
% ONLY WORKS FOR THREE PHASES!!
% 
% INPUTS (all scalars)
% np        sqrt of size of output f
% ffixi     index of phase fraction to fix [optional]
% ffix      value to fix that phase fraction to [optional]
% 
% OUTPUT
% f         output phase fractions [3 x np^2]
% 
% YQW, 25 Jan 2022


if nargin==0, np = 200; end

if nargin<2
    % evenly distributed across three-phase space
    f1  =  linspace(0,1,np);
    f2  =  linspace(0,1,np);

    [f2,f1]  =  meshgrid(f2,f1);

    f2  =  f2(:);
    f1  =  f1(:);
    f3  =  1-f2-f1;

    f   =  [f1,f2,f3].';
    f(:,f(3,:)<0) = nan;

else
    % fix one of the phases, so that we only vary two at a time
    fvary = setdiff(1:3, ffixi);

    f1 = linspace(0, 1-ffix, np^2);

    f          = zeros(3, np^2);
    f(ffixi,:) = ffix.*ones(1,np^2);  % phases with fixed values
    f(fvary,:) = [f1; 1-f1-ffix];   % phases that vary
    f(:,1-f1-ffix<0) = nan;
end


end