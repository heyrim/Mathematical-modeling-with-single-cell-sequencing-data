%% Compute Differentiation Matrix 
function [Dx, Dxx, dx] = FDM_1D( N, dx, norder ) 
    
    N1 = N(1)+1;
    if( nargin < 2 ) 
        dx = 1/N(1); 
    end
    if( nargin < 3 ) 
        norder = 4; 
    end 
    
%%%% 1D 1st deriv 
    Dx = (zeros(N1));   
    if( norder == 4 ) %%%% 4th order 
        Dx(1,1:5) = [-25, 48, -36, 16, -3]./(12*dx); 
        Dx(2,1:5) = [-3,-10, 18, -6, 1]./(12*dx);  
        for n = 3:(N1-2) 
            Dx(n, (n-2):(n+2)) = [1,-8,0,8,-1]./(12*dx); 
        end 
        Dx(N1-1,(N1-4):N1) = [-1, 6, -18, 10, 3]./(12*dx);
        Dx(N1,(N1-4):N1) = [3, -16, 36, -48, 25]./(12*dx); 
        
    elseif( norder == 2 ) %%%% 2nd order          
        Dx(1,1:3) = [-3,4,-1]./(2*dx); 
        for n = 2:(N1-1) 
            Dx(n, (n-1):(n+1)) = [-1,0,1]./(2*dx); 
        end 
        Dx(N1,(N1-2):N1) = [1,-4,3]./(2*dx); 
    end 
    
    Dx = sparse( Dx );      
   
%%%% 1D 2nd deriv - 4th order 
    Dxx = (zeros(N1));   
    if( norder == 4 ) %%%% 4th order     
        Dxx(1,1:5) = [35, -104, 114, -56, 11]./(12*dx^2); 
        Dxx(2,1:5) =   [11, -20, 6, 4, -1]./(12*dx^2); 
        for n = 3:(N1-2) 
            Dxx(n, (n-2):(n+2)) = [-1, 16, -30, 16, -1]./(12*dx^2); 
        end 
        Dxx(N1-1,(N1-4):N1) =  [-1, 4, 6, -20,11]./(12*dx^2); 
        Dxx(N1,(N1-4):N1) = [11, -56, 114, -104, 35]./(12*dx^2); 
    
    elseif( norder == 2 ) %%%% 2nd order          
        Dxx(1,1:4) = [2	-5 4 -1]./(2*dx); 
%         Dxx(2,2:5) = ; 
        for n = 2:(N1-1) 
            Dxx(n, (n-1):(n+1)) = [1,-2,1]./(dx^2); 
        end 
%         Dxx(N1-1,(N1-3):N1) = [-1 4 -5 2]./(2*dx); 
        Dxx(N1,(N1-3):N1) = [-1 4 -5 2]./(2*dx); 

    end    
    
    Dxx = sparse( Dxx ); 
        
    
end 
    
