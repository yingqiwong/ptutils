function [t, xd, zd] = CalcBoxDrift (folder, RunID, t, x, z, f, ustar, wstar)
%
% [t, xd, zd] = RmBoxDriftFromPhaseFrac (folder, RunID)
% 
% calculates the box drift
% 
% YQW, 21 October 2021


fp = GetOutputMatFiles(folder, RunID);
load(fp, 'D');

if nargin<4
    [t,x,z,~,ustar,wstar] = ExtractFieldwTime(folder, RunID, {'f','ustar','wstar'});
end

% calculate mean velocity of reference field
ud   = squeeze(mean(ustar,[2,3]));
wd   = squeeze(mean(wstar,[2,3]));

% calculate the drift in x and z
dx = cumsum([0;ud(2:end).*diff(t)]);
dz = cumsum([0;wd(2:end).*diff(t)]);

% multiple of domain sizes moved
ndx = fix(dx./D(2));
ndz = fix(dz./D(1));

% adjust positions by drift but add back ndz*D if moved by more than 1 domain length
% so that the loop below checks the correct limits.
xd = x - dx + ndx*D(2);
zd = z - dz + ndz*D(1);

end