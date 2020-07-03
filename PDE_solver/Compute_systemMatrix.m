%% Compute PDE Grid and System matrix 
function Compute_systemMatrix( Nin ) 
global N N1 Dxn Dxxn Dxm Dxxm dx xx yy 

N = Nin; 
xx = [0:1/N(1):1]';  yy = [0:1/N(2):1]'; 
dx = [1/N(1), 1/N(2)];
N1 =  N+1; 
    
nN = length( N ); 
Dx = cell(nN,1);  Dxx = cell(nN,1); 
for n = 1:nN 
    [Dx{n}, Dxx{n}] = FDM_1D( N(n)); 
end 
Dxn  = Dx{1};   Dxm  = Dx{2}; 
Dxxn = Dxx{1};  Dxxm = Dxx{2};  

end 
