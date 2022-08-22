function [t, xdc, zdc, f, ustar, wstar, varargout] = RmBoxDrift (folder, RunID, t, x, z, f, ustar, wstar, varmat)
%
% [t, xdc, zdc, fc] = RmBoxDriftFromPhaseFrac (folder, RunID)
% 
% example
% [t, xdc, zdc, fc] = RmBoxDriftFromPhaseFrac('../out/', 'olv10_plg10_bas80');
% 
% calculates the box drift and removes it from the phase fraction field
%
% INPUTS
% folder  	folder name where the simulation is stored
% RunID 	name of simulation
% varnames 	cell vector of variable names [Nvar x 1]
% ti        indices to extract
%
% OUTPUTS
% t         time of simulation [Nf x 1]
% x         x positions of simulation [Nx x 1]
% varmat    cell of variables {Nvar x 1}
%
% YQW, 21 October 2021


fp = GetOutputMatFiles(folder, RunID);
load(fp, 'D');

if nargin<4
    [t,x,z,f,ustar,wstar] = ExtractFieldwTime(folder, RunID, {'f','ustar','wstar'});
end

[t, xd, zd] = CalcBoxDrift(folder, RunID, t, x, z, f, ustar, wstar);

if nargin<9, varmat = {}; end

% now correct position vectors and rearrange matrices so that they lie within [-D/2, D/2].
% initialise
Nt  = length(t);
xdc = xd; 
zdc = zd; 
v0  = [f; ustar; wstar; varmat(:)];  v = v0;
Nv  = length(v);

for ti = 2:Nt
    
    % account for vert and horizontal box drift. Note that you can't have
    % both below & above be not empty, or both left & right be not empty.
    % The box would then not be drifting at one constant speed
    
    % vertical box drift
    below = find(zd(ti,:)<-D(1)/2);
    above = find(zd(ti,:)> D(1)/2,1);
    
    if isempty(below) && isempty(above)
        % limited box drift (within one grid length)
        
    elseif ~isempty(below) && isempty(above)
        % drifting downwards
        zdc(ti,:)    = [zd(ti,below(end)+1:end), zd(ti,below)+D(1)];
        for vi = 1:Nv
            v{vi}(:,:,:,ti,:) = [v0{vi}(:,below(end)+1:end,:,ti,:), v0{vi}(:,below,:,ti,:)];
        end
        
    elseif ~isempty(above) && isempty(below)
        % drifting upwards
        zdc(ti,:)    = [zd(ti,above:end)-D(1) ,zd(ti,1:(above-1))];
        for vi = 1:Nv
            v{vi}(:,:,:,ti,:) = [v0{vi}(:,above:end,:,ti,:), v0{vi}(:,1:(above-1),:,ti,:)];
        end
        
    else
        error('!!!something wrong with vertical box drift!!!');
    end
        
    
    
    % horizontal box drift
    left  = find(xd(ti,:)<-D(2)/2);
    right = find(xd(ti,:)> D(2)/2,1);
    
    if isempty(left) && isempty(right)
        % limited box drift (within one grid length)
    
    elseif ~isempty(left) && isempty(right)
        % drifting leftwards
        xdc(ti,:)    = [xd(ti,left(end)+1:end), xd(ti,left)+D(2)];
        for vi = 1:Nv
            v{vi}(:,:,:,ti,:) = cat(3, v0{vi}(:,:,left(end)+1:end,ti,:), v0{vi}(:,:,left,ti,:));
        end
        
    elseif ~isempty(right) && isempty(left)
        % drifting rightwards
        xdc(ti,:)    = [xd(ti,right:end)-D(2), xd(ti,1:(right-1))];
        for vi = 1:Nv
            v{vi}(:,:,:,ti,:) = cat(3, v0{vi}(:,:,right:end,ti,:), v0{vi}(:,:,1:(right-1),ti,:));
        end
    else
        error('!!!something wrong with horizontal box drift!!!');
    end
end

f = v{1}; ustar = v{2}; wstar = v{3};
v(1:3) = [];

if nargout==7 && length(v)>1
    varargout = {v};
else
    varargout = v;
end

end



