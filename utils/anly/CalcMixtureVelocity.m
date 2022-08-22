function [uvol, wvol] = CalcMixtureVelocity (folder, RunID, f, u, w)
%
% [uvol, wvol] = CalcMixtureVelocity (folder, RunID, f, u, w)
% 
% calculates the volumetric bulk velocity fields for the mixture,
% i.e. vvol = sum_i phi^i v^i
% 
% YQW, 10 Dec 2021



if nargin<3
    [~,~,~,f,u,w] = ExtractFieldwTime(folder, RunID, {'f','u','w'});
end

fvx  = GetCellFaceVals(f,3,'periodic');
fvz  = GetCellFaceVals(f,2,'periodic');

uvol = sum(fvx.*u,1);
wvol = sum(fvz.*w,1);


end