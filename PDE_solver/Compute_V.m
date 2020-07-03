%% Compute Advection to capture homeostasis distribution 
function [V] = Compute_V( uSH )
global cDiff N N1 dx ntest 

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
