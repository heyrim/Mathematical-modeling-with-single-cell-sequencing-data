function Compute_systemCoeff( ntest, uhmsts ) 
global N1 cRct cAdv cDth cA cS V 

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


set_system;


[V] = Compute_V( uhmsts, ntest );


cS = 1-1/(1+exp( -(sum(nCellSave(indSig,1)) - barv)*cK ));


end 


function [V] = Compute_V( uSH, ntest )
global cDiff N N1 dx

ulog = -log( uSH );

Dxn = Compute_Dx2nd( N(1)+1, dx(1) );

if( N(2) == 1 )
    %% 1D
    V = -cDiff * Dxn*ulog;
    
else

    [XX,YY] = meshgrid( 1:N1(2), 1:N1(1) );
    ind = find( ulog ~= inf ) ;
    
    ind = find( uSH > max( uSH(:) )*0.0001 );
    F = scatteredInterpolant(XX(ind), YY(ind), ulog(ind), 'natural', 'none');
    ulog = F(XX,YY); 
    if( ntest == 0 ); ulogMax = 3.5; elseif( ntest == 1 ); ulogMax = 5.5; end 
    ulog(isnan(ulog)) = ulogMax; 
    ulog( ulog>ulogMax ) = ulogMax; 
    
    Dxm = Compute_Dx2nd( N(2)+1, dx(2) );
    
    V{1} = -cDiff * Dxn*ulog;
    V{2} = -cDiff * ulog*Dxm';
    
end

end

