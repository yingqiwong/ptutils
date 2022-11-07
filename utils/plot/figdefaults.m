
function [hAx,pos,ax] = figdefaults (Nrow, Ncol, varargin)
% 
% [hAx,pos] = figdefaults (Nrow, Ncol, varargin)
% 
% function to specify the axis size and margins for figures in paper 


ax.height = 6.00; 
ax.bot    = 1.00; 
ax.top    = 1.00;
ax.gaph   = 1.00; 

ax.totalwidth = 20; % paper width
ax.left   = 2.00;
ax.right  = 0.50;
ax.gapw   = 1.00;

% allow structure alteration
args = reshape(varargin, 2, []);
for ia = 1:size(args,2)
    ax.(args{1,ia}) = args{2,ia};
end

if ~any(strcmp(varargin, 'width'))
    ax.width  = (ax.totalwidth - ax.left - ax.right - (Ncol-1)*ax.gapw)/Ncol;
end

ax.fh = ax.bot  + Nrow*ax.height + (Nrow-1)*ax.gaph + ax.top;
ax.fw = ax.left + Ncol*ax.width  + (Ncol-1)*ax.gapw + ax.right; 
%^^ needed because sometimes I specify width instead of totalwidth

set(gcf,'Units','centimeters','Position',[5 20 ax.fw ax.fh]);
set(gcf,'PaperUnits','Centimeters','PaperPosition',[0 0 ax.fw ax.fh],'PaperSize',[ax.fw ax.fh]);
set(gcf,'defaultaxesfontsize',10);
set(gcf,'defaultlinelinewidth',1);

if Nrow>1 || Ncol>1
    [hAx,pos] = tight_subplot(Nrow,Ncol,[ax.gaph/ax.fh,ax.gapw/ax.fw], [ax.bot,ax.top]/ax.fh, [ax.left,ax.right]/ax.fw);
elseif Nrow==1 && Ncol==1
    hAx = axes('Units','centimeters','Position',[ax.left,ax.bot,ax.width,ax.height]);
end


end
