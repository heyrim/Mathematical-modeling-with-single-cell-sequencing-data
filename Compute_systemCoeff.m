function Compute_systemCoeff( ntestin, uhmsts ) 
global N1 cRct cAdv cDth V2 V indSig cS ntest  

ntest = ntestin; 
if( nargin < 3 ) 
    nAML = 0; 
end 

% Compute Reaction, Self_renewal_fraction, Logistic death, Advection 2  
if( ntest == 0 ) 
    [cRct, cAdv, cDth, V, V2, indSig] = PrePocess_Nestorowa(uhmsts);

elseif( ntest == 1 ) 
    [cRct, cAdv, cDth, V, V2, indSig] = PrePocess_Paul(uhmsts);

end 


%%% Assign Matrix for AML 
Vaml = cell(1,2); 
for n = 1:2; Vaml{n} = zeros( N1(1), N1(2) ); Vaml{n} = sparse( Vaml{n} ); end
addcRct = zeros( N1(1), N1(2) ); addcRct = sparse( addcRct );


set_system_paramters;





end 


