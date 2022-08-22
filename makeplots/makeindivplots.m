

Addpaths

clear all;

set(0, 'defaultlinelinewidth', 2);
set(0, 'defaultaxesfontsize', 18);

set(0, 'defaultAxesTickLabelInterpreter', 'latex'); 
set(0, 'defaultLegendInterpreter'       , 'latex');
set(0, 'defaultTextInterpreter'         , 'latex'); 

beep off

%%
folder    = '../../pantarhei/out/';
RunID = 'olvbas';
fignameprefix = [folder, RunID '/' RunID];

[fp,fn,ft] = GetOutputMatFiles(folder, RunID);

%% load variables

load(fp, 'rho0','h','N','D','grav','delta0','w0');
[t,x,z,f,u,w] = ExtractFieldwTime(folder, RunID, {'f','u','w'});

%%

PlayFieldwTime(folder, RunID, {'f'}, {[]}, 'save',1, 'Nstd', 4);
PlayPhaseFracwTime(folder, RunID);

[t, zp, v] = GetVertProfiles (folder, RunID);
PlotVertProfiles(folder, RunID, {'f'}, [], 'zdsc',true, 'xind', 2);

PlotFieldwTime(folder, RunID, 'f',[],[],'Nstd',5', 'save', 1);

PlotFieldVectors(folder, RunID, {'f','ustar','wstar'}, [], [], [], [], 'save', 1, 'Nplt', 3);

CalcPrincipalStress(folder, RunID, 1, [1,6,11]);

AlternateFrames(folder, RunID, {'f'});















