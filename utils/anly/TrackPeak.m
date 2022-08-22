function [t, fPeak, zPeak, vPeak] = TrackPeak (folder, RunID, gsp, varargin)
%
% [t, fPeak, zPeak, vPeak] = TrackPeak (folder, RunID, gsp, varargin)
%
% example: [t, fPeak, zPeak, vPeak] = TrackPeak('../out/', 'olv20_bas80', [0.25,0.75], varargin)
%
% Track migration of the Gaussian peak from a simulation
% specified by [folder, RunID]
%
% INPUTS
% folder  	folder name where the simulation is stored
% RunID 	name of simulation
% gsp       location of gaussian peak as a fraction of N [Ngsp x 1]
% varargin  options
%
% OUTPUTS
% t         simulation times [Nf x 1]
% fPeak     phase fraction at peak [Nf x Ngsp]
% zPeak     z position of peak as multiple of seg-comp length [Nf x Ngsp]
% vPeak     velocity of peak as multiple of darcy speed [Nf x Ngsp]
%
% YQW, 3 March 2021


%  get output file names
[fp, fn] = GetOutputMatFiles(folder, RunID);
load(fp, 'N', 'D', 'delta0', 'w0', 'BC');

% get lengths
Nf     = length(fn);
Ngsp   = length(gsp);

% plotting options
opt    = defopts(varargin{:});

% initialise
t      = nan(Nf,1   );
fPeak  = nan(Nf,Ngsp);
zPeak  = nan(Nf,Ngsp);

% collect peak info
for fi = 1:Nf
    load(fn{fi}, 'time','f','x','z');
    t(fi) = time;
    
    % collect phase fraction and position at peaks
    for gi = 1:Ngsp
        cut_bc = round(opt.bcind*N);
        [fPeak(fi,gi), zi] = max(f(opt.fi,(cut_bc+1):(end-cut_bc),round(gsp(gi)*N)));
        zPeak(fi,gi)  = z(cut_bc+zi);
    end
end

% account for periodic bcs
if strcmp(BC, 'periodic')
    ziflip = find(zPeak<0 & t>0, 1);
    zPeak(ziflip:end) = z(end) + (zPeak(ziflip:end) - z(1));
end

% normalise z by seg-comp length
zPeak = zPeak./delta0(opt.segind);

% Darcy time scale
tau = delta0(opt.segind)./w0(opt.segind);

% get velocity
vPeak = diff(zPeak)./diff(t./tau);


if opt.plot
    figure;
    set(gcf,'Position',[400,400,1200,300],'defaultlinemarkersize',12);
    tiledlayout(1,3,'TileSpacing','compact');
    
    nexttile;
    plot(t/tau, fPeak, '+:');
    hold on; plot(xlim, fPeak(1)*ones(1,2), 'k:'); hold off;
    xlabel('time [x $\tau_{Darcy}$]'); ylabel('$\delta f^{bas}$');
    title('(a) Basalt frac diff at peak');
    
    nexttile;
    plot(t/tau, zPeak, '+:');
    xlabel('time [x $\tau_{Darcy}$]'); ylabel('z [x max $\delta_{sc}$]');
    title('(b) Location of peak');
    
    nexttile;
    plot(t(2:end)/tau, vPeak, '+:');
    xlabel('time [x $\tau_{Darcy}$]'); ylabel('$v$ [x $w_{Darcy}$]');
    title('(c) Velocity of peak');

    if opt.save
        SaveFigure([folder RunID,'/',RunID,'_gpeakmigration']);
    end
end
end

function opt = defopts (varargin)

opt.fi    = 2;       % phase where gaussian is a peak
opt.bcind = 0;       % fraction of bc to cut off

opt.segind = 3;         % linear index of correct delta0, w0

opt.plot = true;    % plot?
opt.save = false;    % save plot?

% allow structure alteration
args = reshape(varargin, 2, []);
for ia = 1:size(args,2)
    opt.(args{1,ia}) = args{2,ia};
end

end