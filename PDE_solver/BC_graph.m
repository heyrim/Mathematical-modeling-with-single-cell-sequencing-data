function  U = BC_graph( U, time ) 
global BC 
   
%% Dirichlet 
%     U(1,1) = BC.Uinit(1) ; 

%% Neunmann 
    U(1,BC.noin) = U(2,BC.noin); 
%     U(1,BC.Nnmn) = U(2,BC.Nnmn); 
    
    %%%% boundary with no out flow 
    for nn = BC.noout 
      ind =  BC.indin(nn,1:BC.Nindin(nn));   %find( edge(:,2) == nn ); 
%       if( sum(BC.ratein(ind)) )
        tmp = sum( (48*U(end-1,ind)-36*U(end-2,ind)+16*U(end-3,ind)-3*U(end-4,ind)) .*BC.ratein(ind)' , 2 ); 
        tmp = tmp / 25; %/sum(BC.ratein(ind)); 
%       else;   tmp = 0;       
%       end 
      U(end,ind) = tmp; 
%       U(end,ind) = min( tmp, BC.Diclt(nn)); 
%         U(end,ind) = (4*U(end-1,ind)-U(end-2,ind))/3 ;     
    %%%%% although no out.... to other no outs 
      ind =  BC.indout(nn,1:BC.Nindout(nn));   %find( edge(:,2) == nn ); 
      U(1,ind) = tmp; 
    end 
    
   %%%% boundary with in and out flow 
   for nn = BC.inNout 
      indout = BC.indout(nn,1:BC.Nindout(nn));        
      indin  = BC.indin( nn,1:BC.Nindin(nn)); 
      
      tmp = sum( (48*U(2,indout)-36*U(3,indout)+16*U(4,indout)-3*U(5,indout)) .*BC.rateout(indout)' , 2 ) ;  
      tmp = tmp - sum( (-48*U(end-1,indin)+36*U(end-2,indin)-16*U(end-3,indin)+3*U(end-4,indin)) .*BC.ratein(indin)', 2 ) ; 
      tmp = tmp / 50; 
      
      % Diff 
%       tmp = sum( (48*U(2,indout)-36*U(3,indout)+16*U(4,indout)-3*U(5,indout)) .*cDiff(indout)' , 2 ) ;  
%       tmp = tmp - sum( (-48*U(end-1,indin)+36*U(end-2,indin)-16*U(end-3,indin)+3*U(end-4,indin)) .*cDiff(indin)', 2 ) ; 
%       tmp = tmp / 25 / (sum(cDiff(indout))+sum(cDiff(indin))) ; 
      tmp(tmp<0)=0; 
%       disp( abs( sum(cDiff(indout))-sum(cDiff(indin) ) ) ); 
            
      U(1,indout) = tmp; U(end,indin) = tmp; 
    end 
    
end 