function [t, xd, zd] = CalcBoxDrift (folder, RunID, t, x, f, u, w)
%
% [t, xd, zd] = RmBoxDriftFromPhaseFrac (folder, RunID)
% 
% calculates the box drift
% 
% YQW, 21 October 2021


fp = GetOutputMatFiles(folder, RunID);
load(fp, 'D');

if nargin<3
    [t,x,f,u,w] = ExtractFieldwTime(folder, RunID, {'f','u','w'});
end

fvx  = GetCellFaceVals(f,3,'periodic');
fvz  = GetCellFaceVals(f,2,'periodic');

uvol = sum(fvx.*u,1);
wvol = sum(fvz.*w,1);

% calculate mean velocity of reference field
ud   = squeeze(mean(uvol,[2,3]));
wd   = squeeze(mean(wvol,[2,3]));

% calculate the drift in x and z
dx = cumsum([0;ud(2:end).*diff(t)]);
dz = cumsum([0;wd(2:end).*diff(t)]);

% multiple of domain sizes moved
ndx = fix(dx./D);
ndz = fix(dz./D);

% adjust positions by drift but add back ndz*D if moved by more than 1 domain length
% so that the loop below checks the correct limits.
xd = x - dx + ndx*D;
zd = x - dz + ndz*D;

end