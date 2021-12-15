function [ScMod, betaMat] = PlotSuccessfulModels (folder, RunIDVec, Nout)
% 
% [ScMod, betaMat] = PlotSuccessfulModels (folder, RunIDVec, Nout)
% 
% plots whether the simulations of specified phase fractions were
% successful (i.e. ran to specified number of time steps).
% Right now plotting only works for 2, 3 phases.
% 
% INPUTS
% folder  	folder name where the simulations are stored
% RunIDVec 	vector of simulation names [Nruns x 1]
% Nout      number of time steps you expect [1, optional]


if nargin<3, Nout = 1; end


[PHS, NPHS] = ExtractPhsNames(folder, RunIDVec{1});
 Nruns      = length(RunIDVec);

f0Mat   = nan  (Nruns, NPHS);
betaMat = nan  (Nruns, 1   );
ScMod   = false(Nruns, 1   );

for ri = 1:length(RunIDVec)
    
    [fn, fp, ~] = GetOutputMatFiles(folder, RunIDVec{ri});
    
    % load params file
    load(fp, 'f0', 'beta');
    f0Mat(ri,:) = f0;
    betaMat(ri) = beta;
    
    if length(fn) >= Nout
        ScMod(ri) = true;
    end
end




% plot results
if NPHS == 2
    figure;
    lspec = {'bo','MarkerSize',10};
    plot(f0Mat( ScMod,1), betaMat( ScMod), lspec{:}, 'MarkerFaceColor', 'b');
    hold on;
    plot(f0Mat(~ScMod,1), betaMat(~ScMod), lspec{:});
    hold off;
    xlabel(PHS{1}); ylabel(PHS{2});
    
elseif NPHS == 3
    
    figure;
    set(gcf,'Position',[500,500,600,350]);
    [X,Y] = PlotTernaryBG(f0Mat, PHS);
    hold on;
    scatter(X( ScMod), Y( ScMod), 200, betaMat(ScMod), 'filled');
    scatter(X(~ScMod), Y(~ScMod), 200, 'k', 'LineWidth', 2);
    hold off;
    caxis([0, max(betaMat)]);
    cb = colorbar;     
    cb.Location = 'EastOutside';
    cb.Position(1) = cb.Position(1) + 0.1;
    cb.TickLabelInterpreter = 'latex';
    cb.Ticks = min(betaMat):0.1:max(betaMat);
    ht = title(cb, '\beta', 'FontSize', 20);     
end
end