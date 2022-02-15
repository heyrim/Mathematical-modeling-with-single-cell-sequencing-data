function [U] = IC_graph( ninit, usol, ntest ) 
global N1 nEdge nVrtx edge 
 
%     if( nVrtx < 10 ) 
%     load( '180920_Nestorowa_G8.mat' ); 
%     edge = edge( [1,2,3,5,7,9,11], : );
 
    nEdge = size( edge, 1 ); 
    nVrtx = max( max( edge ) ); 
    U = zeros( N1, nEdge ); 
    
    xx = [0:1/(N1(1)-1):1]; 
    sig = 0.3; 
    
%%%% Stem cell concentrated initial condition 
%%%% Approximate fraction on each node with 1D normal distributions   
    if( ntest == 0 ) % Nestorowa 
        
        %%%% IC cluster fraction from continuum initial condition  
        nIC = [1.2087e-01   3.6031e-01   7.2293e-03   1.5514e-05   9.5594e-04   8.6548e-04   3.9457e-03   5.8079e-03]; 

        U(:,1) = Normal( 0, sig, xx )*nIC(2)*2;
        ind = find( edge(:,1) == edge(1,2) ); 
        for n = ind' 
          U(:,n) = Normal( 0, sig, [xx(end)+xx] )*nIC(2)*2; 
        end 
        for n = setdiff( [2:length(edge)], ind ) 
            U(:,n) = U(end,ind(1)); 
        end 

    elseif( ntest==1 ) % Paul 
        xtmp = Normal( 0, sig, xx ); xtmp2 = Normal( 0, sig, [xx(end)+xx] ); 
        ind = find( edge(:,1)==6 )' ; 
        for n = ind 
            U(:,n) = xtmp; 
        end 
        for nn = ind 
          for n = find(edge(:,1)==edge(nn,2))' 
            U(:,n) = xtmp2; 
          end 
        end 
        U(:,4) = xtmp(end); % [6,4] long
        U(:,7) = [xtmp(1:2:end),xtmp2(3:2:end)]; % [1,7] no flow in the middle 
        for n = find( U(1,:) == 0 )'
            U(:,n) = xtmp2(end);  
        end 
        
    end 

    [~, nTotCell] = Compute_nCell_graph( U ); 
    U = U / nTotCell *ninit; 

%%%% Add homeostasis distribution 
    if( ~isempty(usol) ) 
        U = U + usol; 
    end 
    
end 
