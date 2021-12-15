function [KvMat, KfMat, CvMat, CfMat, dscMat, XfMat, thtvMat, thtfMat] = CalcPermFromSims (folder, RunID)
%
% [XfMat, thtvMat, thtfMat, KvMat, KfMat, CvMat, CfMat] = CalcPermissions (folder, RunID, inds)
%
% this function collects connectivity, permissions and coefficients from
% a simulation specified by [folder, RunID]
%
% INPUTS
% folder            folder name where the simulation is stored
% RunID             name of simulation
%
% OUTPUTS
% KvMat  , KfMat    momentum and volume flux coeffs [NPHS x N x N x Nf]
% CvMat  , CfMat    momentum and volume transfer coeffs [NPHS x N x N x Nf]
% dscMat            segregation-compaction lengths [NPHS x NPHS x N x N x Nf]
% XfMat             connectivity matrix [NPHS x NPHS x N x N x Nf]
% thtvMat, thtfMat  momentum and volume permissions [NPHS x N x N x Nf]
%
%
% YQW, 9 Mar 2021


% extract file names
[fn, fp] = GetOutputMatFiles(folder, RunID);
load(fp, 'N','NPHS','A','B','C','d0','eta0','thtlim');

% get diffusivity contrasts
kv = eta0;          % momentum diffusivity
kf = d0.^2./eta0;   % volume diffusivity
Mv = kv.'./kv;      % momentum diffusivity ratios
Mf = kf.'./kf;      % volume diffusivity ratios

% number of time steps
Nf = length(fn);

% initialize matrices
XfMat   = zeros(NPHS, NPHS, N, N, Nf);  dscMat  = XfMat;
thtvMat = zeros(      NPHS, N, N, Nf);  thtfMat = thtvMat;
KvMat   = thtvMat; KfMat = thtvMat; CvMat = thtvMat; CfMat = thtvMat;

for fi = 1:Nf
    load(fn{fi}, 'f', 'delta', 'Kv', 'Kf', 'Cf', 'Cv');
    
    % assign permission outputs
    KvMat (  :,:,:,fi) = Kv;
    KfMat (  :,:,:,fi) = Kf;
    CvMat (  :,:,:,fi) = Cv;
    CfMat (  :,:,:,fi) = Cf;
    dscMat(:,:,:,:,fi) = delta;
    
    if nargout>5
        % get permission weights
        F  = permute(repmat(f,1,1,1,NPHS),[4,1,2,3]);
        Sf = (F./B).^(1./C);  Sf = Sf./sum(Sf,2);
        Xf = sum(A.*Sf,2).*F + (1-sum(A.*Sf,2)).*Sf;
        
        % get momentum and volume permissions
        thtv = squeeze(prod(Mv.^Xf,2));  gmv = geomean(geomean(thtv,3),2);
        thtf = squeeze(prod(Mf.^Xf,2));  gmf = geomean(geomean(thtf,3),2);
        thtv = 1./(1./thtv + 1./(thtlim.^0.5*gmv)) + (gmv/thtlim.^0.5);
        thtf = 1./(1./thtf + 1./(thtlim.^0.5*gmf)) + (gmf/thtlim.^0.5);
        
        % assign outputs
        XfMat  (:,:,:,:,fi) = Xf;
        thtvMat(  :,:,:,fi) = thtv;
        thtfMat(  :,:,:,fi) = thtf;
    end
    
end


end
