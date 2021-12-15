function [t, xdc, zdc, f, u, w, varargout] = RmBoxDrift (folder, RunID, t, x, f, u, w, varmat)
%
% [xdc, zdc, fc] = RmBoxDriftFromPhaseFrac (folder, RunID)
% 
% calculates the box drift and removes it from the phase fraction field
% 
% YQW, 21 October 2021


fp = GetOutputMatFiles(folder, RunID);
load(fp, 'D');

if nargin<3
    [t,x,f,u,w] = ExtractFieldwTime(folder, RunID, {'f','u','w'});
end

[t, xd, zd] = CalcBoxDrift(folder, RunID, t, x, f, u, w);

if nargin<8, varmat = {}; end

% now correct position vectors and rearrange matrices so that they lie within [-D/2, D/2].
% initialise
Nt  = length(t);
xdc = xd; 
zdc = zd; 
v0  = [f; u; w; varmat(:)];  v = v0;
Nv  = length(v);

for ti = 2:Nt
    
    % account for vert and horizontal box drift. Note that you can't have
    % both below & above be not empty, or both left & right be not empty.
    % The box would then not be drifting at one constant speed
    
    % vertical box drift
    below = find(zd(ti,:)<-D/2);
    above = find(zd(ti,:)> D/2,1);
    
    if isempty(below) && isempty(above)
        % limited box drift (within one grid length)
        
    elseif ~isempty(below) && isempty(above)
        % drifting downwards
        zdc(ti,:)    = [zd(ti,below(end)+1:end), zd(ti,below)+D];
        for vi = 1:Nv
            v{vi}(:,:,:,ti,:) = [v0{vi}(:,below(end)+1:end,:,ti,:), v0{vi}(:,below,:,ti,:)];
        end
        
    elseif ~isempty(above) && isempty(below)
        % drifting upwards
        zdc(ti,:)    = [zd(ti,above:end)-D ,zd(ti,1:(above-1))];
        for vi = 1:Nv
            v{vi}(:,:,:,ti,:) = [v0{vi}(:,above:end,:,ti,:), v0{vi}(:,1:(above-1),:,ti,:)];
        end
        
    else
        error('!!!something wrong with vertical box drift!!!');
    end
        
    
    
    % horizontal box drift
    left  = find(xd(ti,:)<-D/2);
    right = find(xd(ti,:)> D/2,1);
    
    if isempty(left) && isempty(right)
        % limited box drift (within one grid length)
    
    elseif ~isempty(left) && isempty(right)
        % drifting leftwards
        xdc(ti,:)    = [xd(ti,left(end)+1:end), xd(ti,left)+D];
        for vi = 1:Nv
            v{vi}(:,:,:,ti,:) = cat(3, v0{vi}(:,:,left(end)+1:end,ti,:), v0{vi}(:,:,left,ti,:));
        end
        
    elseif ~isempty(right) && isempty(left)
        % drifting rightwards
        xdc(ti,:)    = [xd(ti,right:end)-D ,xd(ti,1:(right-1))];
        for vi = 1:Nv
            v{vi}(:,:,:,ti,:) = cat(3, v0{vi}(:,:,right:end,ti,:), v0{vi}(:,:,1:(right-1),ti,:));
        end
    else
        error('!!!something wrong with horizontal box drift!!!');
    end
end

f = v{1}; u = v{2}; w = v{3};
v(1:3) = [];

if nargout==7 && length(v)>1
    varargout = {v};
else
    varargout = v;
end

end



