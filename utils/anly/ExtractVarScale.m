function [varscale, varmat] = ExtractVarScale (folder, RunIDVec, varnames, val, ti)
%
% [varscale, varmat] = ExtractVarScale(folder, RunIDVec, vars, val)
% 
% example: wscls = ExtractVarScale(folder, RunIDVec, {'w','wsegr','wstar'}, 'mid')
%
% this function extracts a scale for variables from a set of simulations
%
% INPUTS
% folder  	folder name where the simulations are stored
% RunIDVec 	vector of simulation names [Nruns x 1]
% varnames 	cell vector of variable names [Nvar x 1]
% val       what kind of scale to extract: 
%               'max': maximum value
%               'mid': midpoint value
% ti        time index to extract from
%
% OUTPUTS
% varscale  cell vector variable scales {1 x Nvar}
% varmat    extracted variable fields {Nruns x Nvar}
%
% YQW, 7 May 2021
%

if nargin<5, ti = 1; end

% get lengths
Nruns = length(RunIDVec);
Nvars = length(varnames);

%  initialise
varmat   = cell(Nruns, Nvars);
varscale = cell(1    , Nvars);

for ri = 1:Nruns
    
    % collect the variable fields from the file
    [~, x, varmat(ri,:)] = ExtractFieldwTime(folder, RunIDVec{ri}, varnames);
    
    % extract the scale for each variable
    for vi = 1:Nvars
        switch val
            case 'min' % get minimum value
                varscale{vi}(:,ri) = min((varmat{ri,vi}(:,:,:,ti)), [], [2,3]);
            
            case 'max'      % get maximum value
                varscale{vi}(:,ri) = max((varmat{ri,vi}(:,:,:,ti)), [], [2,3]);
                
            case 'mid'      % get value at midpoint
                % get index of midpoint
                Nmid = round(0.5*length(x));
                
                varscale{vi}(:,ri) = varmat{ri,vi}(:,Nmid,Nmid,ti);
        end
    end
end

end
