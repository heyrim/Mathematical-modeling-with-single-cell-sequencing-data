function nC = Compute_cluster( nCell )
global Vw dx ntest 

if( ntest == 0 )
    load( 'data/Nestorowa_weights.mat' )
elseif( ntest == 1 )
    load( 'data/Paul_weight.mat' );  
    
end

nVrtx = size( Vw, 3 );
nC = zeros( nVrtx, 1 );
for n = 1:nVrtx
    nC(n) = sum(sum(nCell .* squeeze( Vw(:,:,n) )))*prod(dx); %dxn*dxm;
end


end