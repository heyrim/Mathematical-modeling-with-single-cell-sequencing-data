function [nCell, nTotCell] = Compute_nCell_graph( U ) 
global N1 edge nEdge int_dx nVrtx wUn 
    Nf = (N1+1)/2; 
    
    ind_dxh = int_dx(1:Nf); ind_dxh(end) = ind_dxh(end)/2; 
    nCell = zeros( nVrtx, 1 ); 
    for n = 1:nEdge 
        wtmp = linspace( wUn(1,n), wUn(2,n), N1 )'; 
        nCell(edge(n,1)) = nCell(edge(n,1)) + sum(U(1:Nf, n).*ind_dxh.*wtmp(1:Nf)) ; 
        nCell(edge(n,2)) = nCell(edge(n,2)) + sum(U(Nf:end, n).*ind_dxh.*wtmp(Nf:end)) ; 
    end    

    nTotCell = 0; 
    for n = 1:nEdge 
        wtmp = linspace( wUn(1,n), wUn(2,n), N1 )'; 
        nTotCell = nTotCell + int_dx'*(U(:,n).*wtmp);
    end    
    
    
end 