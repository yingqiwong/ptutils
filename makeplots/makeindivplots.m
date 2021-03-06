

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
RunID = 'porousblayer';
fignameprefix = [folder, RunID '/' RunID];

[fp,fn] = GetOutputMatFiles(folder, RunID);

%% load variables

load(fp, 'rho0','h','N','D','grav','delta0','w0');
[t,x,f,u,w] = ExtractFieldwTime(folder, RunID, {'f','u','w'});


%%

PlayFieldwTime(folder, RunID, {'f'}, {[]}, 'save',1, 'Nstd', 4);

PlotVertProfiles(folder, RunID, {'f'}, [], 'xdsc',1,  'save',1);

PlotFieldwTime(folder, RunID, 'f',[],'Nstd',5,'save',1);

PlotFieldVectors(folder, RunID, {'f','ustar','wstar'}, [], [], [], 'save', 1, 'Nplt', 3);
















