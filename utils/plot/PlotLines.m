function [h] = PlotLines (dir, pos, varargin)
% 
% [h] = PlotLines (dir, pos, varargin)
% 
% plots either horizontal or vertical lines across axis limits at positions
% defined by pos
% 
% INPUTS
% dir       direction of lines ['horiz', 'vert']
% pos       position of lines [Nlines x 1 or 1 x Nlines]
% varargin  line specs
% 
% OUTPUTS
% h         handles to lines
% 
% YQW, 25 Feb 2021


Np = length(pos);
h  = zeros(Np,1);

switch dir
    
    case 'horiz'
        
        hold on;
        for pi = 1:Np
            h(pi) = plot(xlim, pos(pi)*ones(1,2), varargin{:});
        end
        hold off;
        
    case 'vert'
        
        hold on;
        for pi = 1:Np
            h(pi) = plot(pos(pi)*ones(1,2), ylim, varargin{:});
        end
        hold off;
end

uistack(h, 'bottom');

end