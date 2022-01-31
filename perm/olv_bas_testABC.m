% make plots for permission functions and segregation-compaction lengths
% for olv-bas system, testing different values for ABC
% YQW, 19 Jan 2022

clear all;
Addpaths

%% material properties

NPHS = 2;                   % number of phases

PHS  = {'olv','bas'};
rho0 = [ 3200; 2700];       % pure-phase densities
eta0 = [1e+18;1e+02];       % pure-phase viscosities
d0   = [ 5e-3; 5e-3];       % characteristic size of local-scale phase constituents

% set permission weight parameters for coefficient closure model
A = [ 0.6945, 0.1832; 0.5360, 0.1834;];  % permission slopes
B = [ 0.6906, 0.3094; 0.9993, 0.0007;];  % permission step locations
C = [ 0.6889, 0.1750; 0.8154, 1.5642;];  % permission step widths

%% initialize phase fractions

N    = 1e6;
fmid = round(0.5*N);

folv = linspace(0,1,N);
fbas = 1 - folv;
f    = [folv; fbas];

%%
HdrStyle = {'Units','normalized','FontSize',20,'HorizontalAlignment','center'};
LabStyle = {'Units','normalized','FontSize',18,'VerticalAlignment' ,'bottom'};
TtlStyle = {'Units','normalized','FontSize',18};
LegStyle = {'Location','southoutside','FontSize',16,'Box','off', 'NumColumns', 2};

%% calculate permissions given different values of elements in A, B, C

figure;
set(gcf,'Position',[400,400,700,400]);
hAx = tight_subplot(1,2,[0.08,0.03],[0.0,0.1],0.08);
colors = lines(2);
lstyle = {':','-','--'};

% values to test to substitute in various elements in A, B, C
vvec = [0.00001,0.0007,0.4];

for i = 1:3

    Ai = A; Bi = B; Ci = C;

    % assign values in vvec to various elements in A, B or C
    %Ai(1,1) = vvec(i); 
    Bi(2,2) = vvec(i); Bi(2,1) = 1-vvec(i);
    %Ci(2,1) = vvec(i);
    
    [dsc, Kv, Kf, Cv, Cf, Xf] = SegCompLength(f, eta0, d0, Ai, Bi, Ci);
    
    % plot permission weights
    axes(hAx(1));
    hss(i) = plot(fbas, squeeze(Xf(1,1,:)), lstyle{i}, 'Color', colors(1,:)); hold on;
    hll(i) = plot(fbas, squeeze(Xf(2,2,:)), lstyle{i}, 'Color', colors(2,:));
    
    axes(hAx(2));
    hsl(i) = plot(fbas, squeeze(Xf(1,2,:)), lstyle{i}, 'Color', colors(1,:)); hold on;
    hls(i) = plot(fbas, squeeze(Xf(2,1,:)), lstyle{i}, 'Color', colors(2,:));
    
end

vvecstr = repmat(strtrim(string(num2str(vvec'))),2,1);

axes(hAx(1));
xlabel('Liquid fraction, $\phi^\ell$');
ylabel('weights, $X_\phi^{ii}$');
set(gca, 'XTick', 0:0.2:1, 'YTick', 0:0.2:1);
ht=title('(a) Intra-phase weights',TtlStyle{:});
LegText = strcat(repelem({'$ss$, '; '$\ell\ell$, '},3,1), vvecstr);
hleg = legend([hss, hll], LegText, LegStyle{:});
title(hleg, 'value of $B^{2,2}$')

axes(hAx(2));
xlabel('Liquid fraction, $\phi^\ell$');
ylabel('weights, $X_\phi^{ik}$');
set(gca,'YAxisLocation','right');
set(gca, 'XTick', 0:0.2:1, 'YTick', 0:0.2:1);
ht=title('(b) Inter-phase weights',TtlStyle{:});
LegText = strcat(repelem({'$s\ell$, '; '$\ell s$, '},3,1), vvecstr);
hleg = legend([hsl, hls], LegText, LegStyle{:});
title(hleg, 'value of $B^{2,2}$')

