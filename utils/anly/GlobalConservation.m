function [t, cons, Nxz] = GlobalConservation (folder, RunID, plot_opt)

% [t, cons] = GlobalConservation (folder, RunID)
%
% example
% [t, cons] = GlobalConservation('../../pantarhei/out/', 'olv_bas')
%
% Tracks global conservation of phase fractions to check accuracy of
% numerical schemes. A value close to 1 means that the phase is conserved
% well. Each phase should be individually conserved since pantarhei has no
% reactions.
%
% INPUTS
% folder  	folder name where the simulation is stored
% RunID 	name of simulation
%
% OUTPUTS
% t         simulation time [Nt x 1]
% cons      sum of f over domain / initial condition [NPHS x Nt]
%
% YQW, 12 Oct 2022.
%

if nargin<3, plot_opt = false; end

% extract simulation results
[t, ~, ~, f] = ExtractFieldwTime(folder, RunID, {'f'});
Nxz = size(f,2)*size(f,3);

% number of phases
NPHS = size(f,1);

% conservation of transported field, i.e. phase fractions.
cons = squeeze(sum(f,[2,3]));

if plot_opt
    % plot
    figure;
    set(gcf,'Position',[400,400,400,500]);
    tiledlayout(2,1,'Padding','tight','TileSpacing','tight');
    
    nexttile(1);
    plot(t, cons./cons(:,1)-1, '+-');
    hold on; plot(xlim, [0,0], 'k-', 'LineWidth', 0.5); hold off;
    xlabel('time [s]');
    ylabel('$\sum \phi^i (t) \, / \, \sum \phi^i (t=0) - 1$');
    legend(strcat('$\phi^', num2str((1:NPHS)'), '$'),'location','best');
    title('global conservation of $\phi$');
    
    
    nexttile(2);
    plot(t, 1/Nxz*(cons - cons(:,1)), '+-');
    hold on; plot(xlim, [0,0], 'k-', 'LineWidth', 0.5); hold off;
    xlabel('time [s]');
    ylabel('$\left[ \sum \phi^i (t) - \sum \phi^i (t=0)\right]/(N^2)$');
end
end