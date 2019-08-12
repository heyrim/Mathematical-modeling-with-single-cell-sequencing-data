function Compute_systemCoeff( ntest, uhmsts ) 
global N1 

if( ntest == 0 ) 
    
    [cRct, cAdv, cDth, cA] = PrePocess_Nestorowa(uhmsts);
    indSig = [3:6,8];

elseif( ntest == 1 ) 

    [cRct, cAdv, cDth, cA] = PrePocess_Paul(uhmsts);
    indSig = [2,3,4,8,9];

end 



cA(isnan(cA)) = 0;
cAvec{1} = squeeze( cA(:,:,1) );
cAvec{2} = squeeze( cA(:,:,2) );


%%% coeff for AML 
Vaml = cAvec;
for n = 1:2; Vaml{n}(:) = 0; Vaml{n} = sparse( Vaml{n} ); end

addcRct = zeros( N1(1), N1(2) ); addcRct(:) = 0; addcRct = sparse( addcRct );


end 
