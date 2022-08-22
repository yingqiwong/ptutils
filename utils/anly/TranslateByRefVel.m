function [t, z, ft, dzi] = TranslateByRefVel (folder, RunID)
%
% [t, x, f, dxi] = TranslateByRefVel(folder, RunID)
% 
% translates the simulation back by max wstar x t
% assumes periodic BCs
% 

fp = GetOutputMatFiles(folder, RunID);

[t, ~, z, f, wstar] = ExtractFieldwTime(folder, RunID, {'f','wstar'});
Nt = length(t);

wmax = squeeze(max(wstar, [], [2,3]));
dt   = [0; diff(t)];
dz   = cumsum(wmax.*dt);

load(fp, 'h','D','N');

dzi   = floor(dz./h);
Nwrap = floor(dz./D(1));

ft = f;

for ti = 1:Nt
    dzoff = [(dzi(ti)-Nwrap(ti)*N+1):N, 1:(dzi(ti)-Nwrap(ti)*N)];
    ft(:,:,:,ti) = f(:,dzoff,:,ti);
end



end
