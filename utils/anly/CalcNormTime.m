function [tend, dt] = CalcNormTime (t, tau)
% calculate the normalised timescale for different models. 


if iscell(t)
    Nrun = length(t);
    
    % check shape of tau
    if size(tau,1)~=Nrun, tau = tau'; end
    Ntau = size(tau,2);
    
    tend = zeros(Nrun,Ntau); dt = tend;
    for ri = 1:Nrun
        tend(ri,:) = t{ri}(end)./tau(ri,:);
        dt  (ri,:) = (t{ri}(2)-t{ri}(1))./tau(ri,:);
    end
else
    tend = t(end)/tau;
    dt   = (t(2)-t(1))/tau;
end
    


end