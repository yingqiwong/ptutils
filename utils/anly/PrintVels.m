function [] = PrintVels (folder, RunID)


vars = {'f', 'wsegr', 'wstar', 'w', 'p', 'pstar', 'pcmpt'};

load([folder RunID '/' RunID '_par.mat'], 'BC');
load([folder RunID '/' RunID '_0.mat'], vars{:});

xmid   =  floor(0.5*size(f,2));

fprintf('\nRunID\t= %s\n', RunID);

fprintf('f0\t= [ ');
fprintf('%.2f ', f(:,1,1));
fprintf(']\n');

fprintf('wstar\t= [ ');
fprintf('%.2e ', wstar(1,xmid,xmid));
fprintf(']\n');

fprintf('wsegr\t= [ ');
fprintf('%.2e ', wsegr(:,xmid,xmid));
fprintf(']\n');

fprintf('w\t= [ ');
fprintf('%.2e ', w(:,xmid,xmid));
fprintf(']\n');

switch BC
    case {'periodic','open'}
        
        bcfrac = 0;
        exbc   = (floor(bcfrac*size(f,2))+1) : floor((1-bcfrac)*size(f,2));
        
        fprintf('pstar\t= [ ');
        fprintf('%.2e ', max(pstar(1,exbc,exbc),[],[2,3]));
        fprintf(']\n');
        
        fprintf('pcmpt\t= [ ');
        fprintf('%.2e ', max(pcmpt(:,exbc,exbc),[],[2,3]));
        fprintf(']\n');
        
        fprintf('p\t= [ ');
        fprintf('%.2e ', max(p(:,exbc,exbc),[],[2,3]));
        fprintf(']\n\n');
        
    case 'closed'
        
        fprintf('pstar\t= [');
        fprintf('%.2e ', pstar(1,1,xmid));
        fprintf('; ');
        fprintf('%.2e ', pstar(1,end,xmid));        
        fprintf(']\n');
        
        fprintf('pcmpt\t= [');
        fprintf('%.2e ', pcmpt(:,1,xmid));
        fprintf('; ');
        fprintf('%.2e ', pcmpt(:,end,xmid));        
        fprintf(']\n');
        
        fprintf('p\t= [');
        fprintf('%.2e ', p(:,1,xmid));
        fprintf('; ');
        fprintf('%.2e ', p(:,end,xmid));        
        fprintf(']\n');
end

fprintf('\n\n');
   
end

