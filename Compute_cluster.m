
function nC = Compute_cluster( nCell, ntest )
global Vw dx

if( ntest == 0 )
    % load( '190226_Nestorowa_weights.mat' ) % for 101x101
    load( 'data/190227_Nestorowa_weights.mat' )
elseif( ntest == 1 )
    load( 'data/190517_Paul_weight.mat' ); % for 101x101
    % change to 121
%     W(:, 21:121, : ) = Vw(:,:,:);
%     W( 1:26, 1:20, 6) = 1;   W(27:43, 1:20, 7 ) = 1;  W(44:55, 1:20, 5 ) = 1;
%     W(56:end,1:20, 4) = 1/2; W(56:end,1:20, 8 ) = 1/2;
    
    Vw = W;
end

nVrtx = size( Vw, 3 );
nC = zeros( nVrtx, 1 );
for n = 1:nVrtx
    nC(n) = sum(sum(nCell .* squeeze( Vw(:,:,n) )))*prod(dx); %dxn*dxm;
end


end