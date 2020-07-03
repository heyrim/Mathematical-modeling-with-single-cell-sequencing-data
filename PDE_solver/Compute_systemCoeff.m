function Compute_systemCoeff( ntestin, uhmsts, nAML ) 
global N1 cRct cAdv cDth V2 V indSig cS ntest  

ntest = ntestin; 
if( nargin < 3 ) 
    nAML = 0; 
end 

if( ntest == 0 ) 
    [cRct, cAdv, cDth, V2, indSig] = PrePocess_Nestorowa(uhmsts);

elseif( ntest == 1 ) 
    [cRct, cAdv, cDth, V2, indSig] = PrePocess_Paul(uhmsts);

end 


%%% coeff for AML 
Vaml = cell(1,2); 
for n = 1:2; Vaml{n} = zeros( N1(1), N1(2) ); Vaml{n} = sparse( Vaml{n} ); end

addcRct = zeros( N1(1), N1(2) ); addcRct = sparse( addcRct );


set_system;

% nCellSave(:,1) = Compute_cluster( nCell, ntest );
% cS = 1-1/(1+exp( -(sum(nCellSave(indSig,1)) - barv)*cK ));

[V] = Compute_V( uhmsts, ntest );


end 


function [V] = Compute_V( uSH, ntest )
global cDiff N N1 dx

ulog = -log( uSH );

Dxn = FDM_1D( N(1), dx(1), 2 );

if( N(2) == 1 )    %% 1D
    
    V = -cDiff * Dxn*ulog;
    
else   %% 2D 
    
    [XX,YY] = meshgrid( 1:N1(2), 1:N1(1) );
%     ind = find( ulog ~= inf ) ;
    
    ind = find( uSH > max( uSH(:) )*0.0001 );
    F = scatteredInterpolant(XX(ind), YY(ind), ulog(ind), 'natural', 'none');
    ulog = F(XX,YY); 
    if( ntest == 0 ); ulogMax = 3.5; elseif( ntest == 1 ); ulogMax = 5.5; end 
    ulog(isnan(ulog)) = ulogMax; 
    ulog( ulog>ulogMax ) = ulogMax; 
    
    Dxm = FDM_1D( N(2), dx(2), 2 );
    
    V{1} = -cDiff * Dxn*ulog;
    V{2} = -cDiff * ulog*Dxm';
    
end

end

