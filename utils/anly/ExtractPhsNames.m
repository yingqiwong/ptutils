function [PHS, NPHS] = ExtractPhsNames (folder, RunID)
% 
% [PHS, NPHS] = ExtractPhsNames (folder, RunID)
% 
% example: [PHS, NPHS] = ExtractPhsNames('../out/', 'olv20_bas80')
% 
% extract phase names from simulation specified by [folder, RunID]. 
% If PHS does not exist, then the output is just {f1, f2, etc.}
% 
% INPUTS
% folder  	folder name where the simulation is stored
% RunID 	name of simulation
% 
% OUTPUTS
% PHS       phase names [NPHS x 1]
% NPHS      number of phases [1]
% 
% YQW, 19 Jan 2021
% 

% try to load PHS from file. May or may not work depending on the file
load([folder RunID '/' RunID '_par.mat'], 'PHS');

if ~exist('PHS','var')
    % if not stored in the file and needs to be extracted from filename
    load([folder RunID '/' RunID '_par.mat'], 'f0');
    
    NPHS = length(f0);
    PHS  = strcat({'f'}, num2str((1:NPHS)', '%d'));

else
    NPHS = length(PHS);
end



end