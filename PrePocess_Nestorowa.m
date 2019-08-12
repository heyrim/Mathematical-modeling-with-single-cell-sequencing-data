function [cPmat, cAmat, cDth, cAvec] = PrePocess_Nestorowa(uhmsts)  
global N 

load( 'Nestorowa2016_scRNAseqData.mat' ) 
cellind = cluster; cellind( cellind == 9 ) = 2; 


%% Mesh 
[YY,XX] = meshgrid([0:(1/N(2)):1], [0:(1/N(1)):1]);

%%%% add zeros at coarse grids 
len = 0.04; 
[Y,X] = meshgrid([len/2:(len):(1-len/2)], [len/2:(len):(1-len/2)]); 
% [Y,X] = meshgrid([0:(len):1], [0:(len):1]); 


data = dc_normalize(dc, 2);
datawzero = [data; X(:),Y(:)];  
M = size(datawzero,1); 


%% Prolif 
%%%% log(2)/ [46*7 8.8*7 4.3 52/24] 
% cProlif = [0.002152   0.01125   0.16120   0.32]; 
%%%% log(2)/ [46*7 8.8*7 12.25 4.3 52/24 24/24] 
cProlif = [0.002152  0.01125  0.05658  0.1612  0.3199  0.6931]; 
cProlif = cProlif( [2 3 4 6] ); 
cellind2state = [2 1 3 4 4 3 2 3]; 

cPdata = cProlif( cellind2state( cellind ) )'; 

%%%% Global data interp with sin(nx) 
% [cP,ndeg] = data_interp( data, cPdata, 2, 3, 1 );  % nspecial=1 for sin(nx) 
% cPmat = polyfitn( cP, ndeg, [XX(:), YY(:)], 1);

%%%post process for negative 
% ii = find( XX<0.5 ); ii = intersect( ii, find( YY<0.5 ) ); ii = intersect( ii, find( cPmat<cProlif(1) ) ) ; ii = setdiff( ii, find( cPmat==0 ) ); 
% cPmat(ii) = cProlif(1); 
% cPmat = reshape( cPmat, N(1)+1, N(2)+1 ); 


cPdata(M) = 0; 

cPmat = zeros( size(X,1), size(X,2) ); 
for n = 1:size(X,1); for m = 1:size(X,2) 
    dd = sqrt( sum( bsxfun(@minus, datawzero ,[X(n,m),Y(n,m)]).^2, 2 ) ); 
    ii = find( dd < len ); 
    cPmat(n,m) = mean(cPdata(ii,:),1); 
end; end 
cPmat = interp2( X', Y', cPmat', XX', YY', 'linear', 0 )'; 

%%%%% plot results 
% figure; surfo( XX, YY, cPmat ) 


%% differentiation fraction 
cAdvec = [0.77 0.7689 0.7359 0.66 0.154]; 
cAdvec = cAdvec(1:4); 
cAdata = cAdvec( cellind2state( cellind ) )'; 

% [cA,ndeg] = data_interp( data, cAdata, 2, 3 );
% cAmat = polyfitn( cA, ndeg, [XX(:), YY(:)]);
% 
% %%%post process for 0.11<cA<1-alpha 
% cAmat(cAmat<0.11)=0.11; cAmat(cAmat>0.77)=0.77; 
% cAmat = reshape( cAmat, N(1)+1, N(2)+1 ); 

cAdata( (length(cAdata)+1):M ) = 1; 

cAmat = zeros( size(X,1), size(X,2) ); 
for n = 1:size(X,1); for m = 1:size(X,2) 
    dd = sqrt( sum( bsxfun(@minus, datawzero ,[X(n,m),Y(n,m)]).^2, 2 ) ); 
    ii = find( dd < len ); 
    cAmat(n,m) = mean(cAdata(ii,:),1); 
end; end 
cAmat = interp2( X', Y', cAmat', XX', YY', 'linear', 1 )'; 

%%%%% plot results 
% figure; surfo( XX, YY, cAmat ) 


%% direction of differentiation 
%%%%% reference find_DPF_vector.m 

load('/Users/heyrim/Downloads/Software/DPT_inMatlab/dpt/examples/190220_DPT.mat') 
dxdynew(isnan(dxdynew)) = 0; 


%%%% normalize to one 
for n = 1:length(dxdynew) 
    ntmp = sqrt(sum(dxdynew(n,:).^2)); 
    if( ntmp );  dxdynew(n,:) = dxdynew(n,:)/ntmp; end; 
end

dxdy = dxdynew;
dxdy(M, 2) = 0; 

dxdyMat = zeros( size(X,1), size(X,2), 2 ); 
for n = 1:size(X,1); for m = 1:size(X,2) 
    dd = sqrt( sum( bsxfun(@minus, datawzero ,[X(n,m),Y(n,m)]).^2, 2 ) ); 
    ii = find( dd < len ); 
    dxdyMat(n,m,:) = mean(dxdy(ii,:),1); 
end; end 
% figure; quiver(X, Y, dxdyMat(:,:,1), dxdyMat(:,:,2) ); 
for n = 1:2 
cAvec(:,:,n) = interp2( X', Y', dxdyMat(:,:,n)', XX', YY', 'linear', 0 )'; 
end 


uhmsts(uhmsts<=eps)=eps;
cDth = 1./uhmsts;


end 

function [data] = addzero( data, add ) 

data = [data; add]; 

end 

function [newvec] = AMLvec()
dd = figure_extract( 6, 1 );
ind = find( cellind == 4 );
dd.x{1} = dd.x{1}(ind); 
dd.y{1} = dd.y{1}(ind);
vec(:,1) = -dd.x{1}' + dd.x{2}' ;
vec(:,2) = -dd.y{1}' + dd.y{2}' ;

datawzero = [data(ind,1),data(ind,2); X(:),Y(:)];
dxdy = vec; dxdy( size(datawzero,1), 2 ) = 0; 

save( '190406_Nstrw_AML_vec.mat', 'XX', 'YY', 'newvec', 'dd' );

end 



