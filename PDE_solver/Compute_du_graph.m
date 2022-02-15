function [du] = Compute_du_graph( U, time ) 
global  Dxx Dx cDiffM cDiff  cS  cAdvM cRctM cDthM cAdvM2 addcRct  
    %% PDF on graph    
    React = ( 1 - cDthM.*U );   React(React<-1)=-1; 
     
    du = cDiff * (Dxx*U)*cDiffM - Dx* (cAdvM.*U) + cRctM.*React.*U; % + cSp*cRctM.*U ; 
    du = du - cS*Dx *( cAdvM2.*U ); 
    du = du + addcRct.*U; 
    
%     du = cDiff * (Dxx*U)*cDiffM  - Dx *( cAdvM2.*U ) + cRctM.*React.*U; 
        
end