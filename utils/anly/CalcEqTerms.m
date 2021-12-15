function [] = CalcEqTerms (folder, RunID, varargin)
% 
% CalcEqTerms (folder, RunID, varargin)
% calculates individual terms of the z momentum and mass balance equations
% 
% example
% CalcEqTerms (folder, RunID, 'ti', 4);
% 
% 
% INPUTS
% folder    folder name where output folder is stored
% RunID     name of the run so that the total path to a mat file is
%               [folder  RunID / RunID_0.mat]
% varargin  plotting options (see defopts)
% 
% 
% DEFAULT OPTIONS
% opt.plot = 1; 
% opt.save = 0;
% opt.ti   = 1;
% 
% YQW, 15 June 2021



[fp,fn] = GetOutputMatFiles(folder, RunID);

opt = defopts(varargin{:});

load(fp, 'h','BC','N','NPHS','grav','rho0');
load(fn{1}, 'x','f');
rho    = rho0.*ones(NPHS,N,N);
rhomix = mean(mean(sum(f.*rho,1)));

% initialise indexing for boundary condition stencils
if     strcmp(BC,'periodic'); ic = [N,1:N,1]; im = [N,1:N]; ip = [1:N,1];
elseif strcmp(BC,'open') || strcmp(BC,'closed'); ic = [1,1:N,N]; im = [1,1:N]; ip = [1:N,N]; end

vars = {'f','p','pstar','pcmpt','u','ustar','usegr','w','wstar','wsegr',...
    'Kv','Cv','Cf','qvxx','qvxz','qvzz','qfx','qfz','Gvz'};
load(fn{opt.ti}, vars{:});


% phase-wise x momentum----------------------------------------------------
fu  = (f(:,:,im)+f(:,:,ip))./2;
Cvu = (Cv(:,:,im)+Cv(:,:,ip))./2;

qxdev = qvxx - f.*p;

dpstarx = -fu.^2./Cvu.*diff(pstar(:,:,ic),1,3)./h;
dpcmptx = -fu   ./Cvu.*diff(pcmpt(:,:,ic),1,3)./h;
dKvDx   = -fu.^2./Cvu.*(diff(qvxz,1,2)./h + diff(qxdev(:,:,ic),1,3)./h);
Qvx     = -fu.^3./Cvu.*((rho(:,:,im)+rho(:,:,ip))./2-rhomix).*grav(2);

vx = {usegr, dpstarx, dpcmptx, dKvDx, Qvx};
vxnames = {'u_\Delta','\phi^2/Cv \times \nabla P^*', ...
    '\phi/Cv \times \nabla P^{cmpt}',...
    '\phi^2/Cv \times \nabla \cdot (Kv D)','\phi^2/Cv \times Qvx'};

if opt.plot, plotterms(folder, RunID, x, vx, 'vx', vxnames, opt); end



% phase-wise z momentum----------------------------------------------------
fw  = (f(:,im,:)+f(:,ip,:))./2;
Cvw = (Cv(:,im,:)+Cv(:,ip,:))./2;

qzdev = qvzz - f.*p;

dpstarz = -fw.^2./Cvw.*diff(pstar(:,ic,:),1,2)./h;
dpcmptz = -fw   ./Cvw.*diff(pcmpt(:,ic,:),1,2)./h;
dKvDz   = -fw.^2./Cvw.*(diff(qvxz,1,3)./h + diff(qzdev(:,ic,:),1,2)./h);
Qvz     = -fw.^3./Cvw.*((rho(:,im,:)+rho(:,ip,:))./2-rhomix).*grav(1);

vz = {wsegr, dpstarz, dpcmptz, dKvDz, Qvz};
vznames = {'w_\Delta','\phi^2/Cv \times \nabla P^*', ...
    '\phi/Cv \times \nabla P^{cmpt}',...
    '\phi^2/Cv \times \nabla \cdot (Kv D)','\phi^2/Cv \times Qvz'};

if opt.plot, plotterms(folder, RunID, x, vz, 'vz', vznames, opt); end




% phase-wise mass balance--------------------------------------------------
qfxp = qfx - (f(:,:,im)+ f(:,:,ip))./2 .* u;
qfzp = qfz - (f(:,im,:)+ f(:,ip,:))./2 .* w;

dvstar = -f.^2./Cf.*(diff(ustar,1,3)./h + diff(wstar,1,2)./h);
dvsegr = -f   ./Cf.*(diff(usegr,1,3)./h + diff(wsegr,1,2)./h);
dKfdp  = -f.^2./Cf.*(diff( qfxp,1,3)./h + diff( qfzp,1,2)./h);

mass = {pcmpt, dvstar, dvsegr, dKfdp};
massnames = {'p^{cmpt}', '\phi^2/C_\phi \times \nabla \cdot v^*', ...
    '\phi/C_\phi \times \nabla \cdot v_\Delta', ...
    '\phi^2/C_\phi \times \nabla \cdot K_\phi (\nabla P)^{i*}'};

if opt.plot, plotterms(folder, RunID, x, mass, 'mass', massnames, opt); end





% mixture mass balance-----------------------------------------------------
divvstar =     diff(ustar,1,3)./h + diff(wstar,1,2)./h;
divvsegr = sum(diff(usegr,1,3)./h + diff(wsegr,1,2)./h, 1);

mixmass      = {divvstar, divvsegr};
mixmassnames = {'\nabla \cdot v^*', '\sum_i \nabla \cdot v_\delta^i'};

if opt.plot, plotterms(folder, RunID, x, mixmass, 'mixmass', mixmassnames, opt); end



% mixture z-momentum balance-----------------------------------------------
gradpstar   =     diff(pstar(:,:,ic),1,2)./h;
gradpcmpt   = sum(diff(pcmpt(:,:,ic),1,2)./h, 1);
gradpbar    = gradpstar + gradpcmpt;
shearstress = sum( diff(qvxz,1,3)./h + diff(qzdev(:,ic,:),1,2)./h, 1);
Qvz         = (sum(f.*rho,1)-rhomix).*grav(1);

mixmv      = {gradpstar, gradpcmpt, gradpbar, shearstress, Qvz};
mixmvnames = {'\nabla p^*', '\sum_i \nabla p_\Delta^i', '\nabla \bar{p}', ...
    '\sum_i \nabla \cdot K_v^i D^i', '\Delta \bar{\rho} g'};

if opt.plot, plotterms(folder, RunID, x, mixmv, 'mixmv', mixmvnames, opt); end


end


function plotterms (folder, RunID, x, v, eqtype, tnames, opt)

% load colormap
load('../pantarhei/src/ocean.mat', 'ocean');

Nv   = length(v);
NPHS = size(v{1},1);

if nargin<6, tnames = strcat({'t'}, num2str((1:Nv)')); end

% prepare for plotting
TX = {'Interpreter','Latex'}; FS = {'FontSize',18};
TL = {'TickLabelInterpreter','Latex'}; TS = {'FontSize',14};

% initialize figure and axes
fig = figure;
colormap(ocean);
set(fig,'Position',[500,500,300*Nv+40,240*NPHS+60]);
set(fig,'Name',[RunID, '_' eqtype]);

hAx = tight_subplot(NPHS,Nv,0.02,[0.03,0.05],0.05);

% plot panels
for m = 1:Nv
    
    for iphs=1:NPHS
        axes(hAx((iphs-1)*Nv+m));
        
        % plot background
        imagesc(x,x,squeeze(v{m}(iphs,:,:)));  
        
        axis xy equal tight;
        xlim([min(x), max(x)]); ylim([min(x), max(x)]);
        
        cb = colorbar; set(cb,TL{:},TS{:});
        set(gca,TL{:},TS{:}); set(gca,'XTickLabel',[]);
        ax((iphs-1)*Nv+m).YAxis.Exponent = 0;
        
        climits = get(gca,'clim');
        if sign(climits(1))~=sign(climits(2)), set(gca,'clim',max(abs(climits)).*[-1,1]); end
        
        if m>1, set(gca,'YTickLabel',[]); end
        if iphs==1, title(['$' tnames{m} '$'],TX{:},FS{:}); end
        
    end
end

if opt.save
    figname = [folder, RunID '/' RunID '_' eqtype 'eqterms'];
    SaveFigure(figname, fig);
end
end

function [opt] = defopts (varargin)

opt.plot = 1; 
opt.save = 0;
opt.ti   = 1;

% allow structure alteration
args = reshape(varargin, 2, []);
for ia = 1:size(args,2)
    opt.(args{1,ia}) = args{2,ia};
end

end





