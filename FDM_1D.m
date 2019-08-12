%% Return Differentiation Matrix 
function [Dx, Dxx, dx] = FDM_1D( N ) 
    
    N1 = N(1)+1;
    dx = 1/N(1); 
    
%%%% 1D 1st deriv - 4th order 
    Dx = (zeros(N1));   
    Dx(1,1:5) = [-25, 48, -36, 16, -3]./(12*dx); 
    Dx(2,1:5) = [-3,-10, 18, -6, 1]./(12*dx);  
    for n = 3:(N1-2) 
        Dx(n, (n-2):(n+2)) = [1,-8,0,8,-1]./(12*dx); 
    end 
    Dx(N1-1,(N1-4):N1) = [-1, 6, -18, 10, 3]./(12*dx);
    Dx(N1,(N1-4):N1) = [3, -16, 36, -48, 25]./(12*dx); 
    Dx = sparse( Dx );        

%%%% 1D 2nd deriv - 4th order 
    Dxx = (zeros(N1));   
    Dxx(1,1:5) = [35, -104, 114, -56, 11]./(12*dx^2); 
    Dxx(2,1:5) =   [11, -20, 6, 4, -1]./(12*dx^2); 
    for n = 3:(N1-2) 
        Dxx(n, (n-2):(n+2)) = [-1, 16, -30, 16, -1]./(12*dx^2); 
    end 
    Dxx(N1-1,(N1-4):N1) =  [-1, 4, 6, -20,11]./(12*dx^2); 
    Dxx(N1,(N1-4):N1) = [11, -56, 114, -104, 35]./(12*dx^2); 
    Dxx = sparse( Dxx ); 
        
end 
