
function [fp, fn, NPHS] = GetOutputMatFiles (folder, RunID)
% 
% [fn, fp, NPHS] = GetOutputMatFiles (folder, RunID)
% 
% example: [fn, fp, NPHS] = GetOutputMatFiles('../out/', 'olv20_bas80')
% 
% extracts the names of simulation and parameter files from a simulation
% specified by [folder, RunID]
% 
% INPUTS
% folder  	folder name where the simulation is stored
% RunID 	name of simulation
% 
% OUTPUTS
% fn        cell vector containing all the names of simulations
% fp        string of parameter file name
% NPHS      number of phases
% 
% YQW, 19 January 2021
% 

% get output mat files and sort in order
outdir	= [folder RunID '/'];
f    	= dir([outdir '*.mat']);

% separate files into simulation and parameter files
fname = strcat(outdir, {f.name}');
fp    = char(fname( contains(fname, '_par')));
fn    =      fname(~contains(fname, '_par'));

% if you only want the parameters file, quit
if nargout==1, return; end

% collect number at end of simulation file name so we can sort in order
ind = zeros(length(fn),1);
for fi = 1:length(fn)
    tmp     = strsplit(fn{fi},'_');
    tmp2    = strsplit(tmp{end}, '.');
    ind(fi) = str2double(tmp2{1});
end

% sort by simulation index
[~,order] = sort(ind,'ascend');
fn        = fn(order);

% check the number of phases
if nargout>2, load(fp, 'NPHS'); end
end
