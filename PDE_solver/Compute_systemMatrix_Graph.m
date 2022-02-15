function Compute_systemMatrix_Graph( Nin )
global N N1 x dx int_dx Dx Dxx 

    N = Nin; 
    N1 = N+1; 
    L1min = 0; L1max = 1; 
    x   = linspace(L1min, L1max, N1(1))';  dx = (L1max-L1min)/N(1); 
    int_dx = ones(N1(1),1)*dx; int_dx(1) = int_dx(1)/2; int_dx(end) = int_dx(end)/2; 
    
   [Dx, Dxx] = FDM_1D( N(1) ); 

end 