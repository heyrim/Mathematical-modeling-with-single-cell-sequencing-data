
function nC = Compute_cluster( nCell )
global Vw dx ntest 

if( ntest == 0 )
    % load( '190226_Nestorowa_weights.mat' ) % for 101x101
    load( 'data/190227_Nestorowa_weights.mat' )
elseif( ntest == 1 )
    load( 'data/190517_Paul_weight.mat' ); % for 101x101
    
end

nVrtx = size( Vw, 3 );
nC = zeros( nVrtx, 1 );
for n = 1:nVrtx
    nC(n) = sum(sum(nCell .* squeeze( Vw(:,:,n) )))*prod(dx); %dxn*dxm;
end


end