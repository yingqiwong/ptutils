function [t, x, ft, dxi] = TranslateByRefVel (folder, RunID)
%
% [t, x, f, dxi] = TranslateByRefVel(folder, RunID)
% 
% translates the simulation back by max wstar x t
% assumes periodic BCs
% 

fp = GetOutputMatFiles(folder, RunID);

[t, x, f, wstar] = ExtractFieldwTime(folder, RunID, {'f','wstar'});
Nt = length(t);

wmax = squeeze(max(wstar, [], [2,3]));
dt = [0; diff(t)];
dx = cumsum(wmax.*dt);

load(fp, 'h','D','N');

dxi   = floor(dx./h);
Nwrap = floor(dx./D);

ft = f;

for ti = 1:Nt
    dxoff = [(dxi(ti)-Nwrap(ti)*N+1):N, 1:(dxi(ti)-Nwrap(ti)*N)];
    ft(:,:,:,ti) = f(:,dxoff,:,ti);
end



end
