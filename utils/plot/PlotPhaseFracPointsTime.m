function [t, fpt] = PlotPhaseFracPointsTime (folder, RunID, xi, zi)
% PlotPhaseFracPointsTime (folder, RunID)
% 
% plots the phase fraction for some points in the domain through time
% 
% INPUTS
% folder    folder name where runs are stored
% RunID     name of run
% xi, zi    x and z indices of specified points to plot
%           should be same size [Nxi x 1 or transpose]
% 
% OUTPUTS
% t         time of solutions (normalized by shortest timescale) [Nf x 1]
% fpt       phase fraction at those points [NPHS x Nxi x Nf]


[fn, fp] = GetOutputMatFiles(folder, RunID);
Nf = length(fn);

% load necessary information from files
load(fp, 'delta0', 'DeltaRho0', 'w0', 'D', 'Kv', 'grav', 'NPHS','N', 'PHS');
if nargin == 2, [xi, zi] = DefaultInds (N); end
Nxi = length(xi);


% get velocity scales and timescales --------------------------------------
% Stokes scale
wSt   = DeltaRho0.*max(abs(grav)).*D.^2./sum(Kv);
tauSt = D./max(wSt(:));

% Darcy timescale
tauDy = max(delta0(:))./max(w0(:));

% choose shortest timescale
tau = min(tauSt, tauDy);


% extract phase fractions at specified points -----------------------------
% convert [zi,xi] into linear indices
linInd = sub2ind([NPHS,N,N], repelem(1:NPHS,Nxi), repmat(zi,1,NPHS), repmat(xi,1,NPHS));

% initialize output matrices
fpt = zeros(NPHS, Nxi, Nf);
t   = zeros(Nf, 1);

for ti = 1:length(fn)
    load(fn{ti}, 'time', 'f','x');
    t(ti) = time;
    fpt(:,:,ti) = reshape(f(linInd),length(xi),NPHS)';
    
    if ti == 1, ft0 = f; end
end

t = t/tau;


% plot --------------------------------------------------------------------
% load colormap
load ../src/ocean.mat;

% prepare for plotting
TX = {'Interpreter','Latex'}; FS = {'FontSize',18};
TL = {'TickLabelInterpreter','Latex'}; TS = {'FontSize',14};
MS = {'MarkerSize', 10};

figure;
set(gcf,'Position',[400,400,1000,NPHS*250]);
colormap(ocean);
hAx = tight_subplot(NPHS,3, [0.1,0.03], [0.1,0.05]);

for iphs = 1:NPHS
    axes(hAx((iphs-1)*3 + 1));
    imagesc(x,x,squeeze(ft0(iphs,:,:)));
    hold on; plot(x(xi), x(zi), 'w+', MS{:}); hold off
    text(x(xi), x(zi), strcat({'\,\, '}, num2str((1:Nxi)', '%d')), 'Color', 'w', 'FontSize',12);
    axis xy equal tight;
    cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:});
    set(gca,'XTick', []);
    title(['$\phi^{' PHS{iphs} '} (t = 0)$'], FS{:});
    
    axes(hAx((iphs-1)*3 + 2));
    hAx((iphs-1)*3 + 2).Position(1) = hAx((iphs-1)*3 + 2).Position(1) - 0.02;
    imagesc(x,x,squeeze(f(iphs,:,:)));
    hold on; plot(x(xi), x(zi), 'w+', MS{:}); hold off
    text(x(xi), x(zi), strcat({'\,\, '}, num2str((1:Nxi)', '%d')), 'Color', 'w', 'FontSize',12);
    axis xy equal tight;
    cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:});
    set(gca,'XTick', []);
    title(['$\phi^{' PHS{iphs} '} (t =$ end)'], FS{:});
    
    axes(hAx((iphs-1)*3 + 3));
    plot(t, squeeze(fpt(iphs,:,:)), '+:', MS{:});
%     colororder(gca, parula(Nxi+1));
    text(t(end)*ones(Nxi,1), fpt(iphs,:,end), strcat({'\,\, '}, num2str((1:Nxi)', '%d')), TS{:});
    ylabel(['$\phi^{' PHS{iphs} '} (t)$']);
    if iphs == NPHS, xlabel('$t/\tau$'); end

end

end

function [xi, zi] = DefaultInds (N)

xi = repmat(round(0.5*N), 1, 3);
zi = [2, round(0.5*N), N-1];

end