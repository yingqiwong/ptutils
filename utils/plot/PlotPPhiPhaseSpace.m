function [fig] = PlotPPhiPhaseSpace (folder, RunID, gi, FileInd, inds)
% 
% [fig] = PlotPPhiPhaseSpace (folder, RunID, gi, pi, fi);
% 
% Use this to plot the p-phi phase space as in Connolly and Podladchikov
% (1998). p is normalized by the compaction pressure (drho x g x delta_c),
% while phi is normalized by the background level.
% 
% INPUTS
% folder    folder name where runs are stored
% RunID     name of run
% gi        x index location to plot [scalar]
% FileInd   file numbers to plot [Nf x 1]
% i         phase indices for f, p respectively
% 
% YQW, 3 Mar 2021


[fn, fp] = GetOutputMatFiles(folder, RunID);
fi = inds(1);   % index of phase for f
pi = inds(2);   % index of phase for p

if isempty(FileInd), FileInd = 1:length(fn); end

load(fp,'h','delta0','f0','rho0','grav','PHS','NPHS');
if ~exist('PHS', 'var'), PHS = strcat({'f'}, num2str((1:NPHS)', '%d')); end

[Nrow,Ncol] = GetSubplotRowCol(length(FileInd));

fig=figure;
tiledlayout(Nrow,Ncol,'TileSpacing','compact');

for ti = FileInd
    load(fn{ti}, 'f', 'pcmpt');
    
    p0  = abs(diff(rho0(inds)).*grav(1).*max(delta0(:)));
    f2z =     f(fi,:,gi)./f0(fi);
    p2z = pcmpt(pi,:,gi)./p0;
    
    dfdz = gradient(f2z, h);
    dpdz = gradient(p2z, h);
    
    nexttile;
    quiver(f2z, p2z, dfdz, dpdz);
    axis equal
    xlabel(['$\phi^{' PHS{fi} '}/\phi^{' PHS{fi}  '}_0$']);
    ylabel(['$p_\textrm{cmpt}^{' PHS{pi} '}/(\Delta\rho g \delta_{sc})$']);
    title(['ti = ' num2str(ti)]);
end



end