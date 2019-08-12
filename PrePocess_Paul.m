function  [cPmat, cAmat, cDth, cAvec] = PrePocess_Paul(uhmsts) 
global N 

load( 'Paul2015_scRNAseqData.mat' );  
cluster(cluster==10) = 6; cellind = cluster; 
   
%% Mesh 
y = [-0.2:1.2/N(2):1]'; 
[YY,XX] = meshgrid(y,[0:(1/N(1)):1]);

%%%% add zeros at coarse grids 
len = 0.04; 
[Y,X] = meshgrid([(-0.2+len/2):(len):(1-len/2)], [len/2:(len):(1-len/2)]); 
% [Y,X] = meshgrid([0:(len):1], [0:(len):1]); 


data = dc;   %_normalize(dc, 2);
datawzero = [data; X(:),Y(:)];  
M = size(datawzero,1);    

%% Prolif 
cProlif = [0.002152  0.01125  0.05658  0.1612  0.3199  0.6931]; 
cProlif = cProlif( [2 3 4 6] ); 
cellind2state = [ 2 4 4 4 3 1 2 3 2 ]; 
cProlif = cProlif( cellind2state ); 
cPdata = cProlif( cellind2state( cellind ) )'; 
cPdata(M) = 0; %%%% fill zero 

cPmat = zeros( size(X,1), size(X,2) ); 
for n = 1:size(X,1); for m = 1:size(X,2) 
    dd = sqrt( sum( bsxfun(@minus, datawzero ,[X(n,m),Y(n,m)]).^2, 2 ) ); 
    ii = find( dd < len ); 
    cPmat(n,m) = mean(cPdata(ii,:),1); 
end; end 
cPmat = interp2( X', Y', cPmat', XX', YY', 'linear', 0 )'; 


%% differentiation fraction 
cAdvec = [0.77 0.7689 0.7359 0.66 0.154]; 
cAdvec = cAdvec(1:4); 
cAdata = cAdvec( cellind2state( cellind ) )'; 
cAdata( (length(cAdata)+1):M ) = 1; %%%% fill 1 
   
cAmat = zeros( size(X,1), size(X,2) ); 
for n = 1:size(X,1); for m = 1:size(X,2) 
    dd = sqrt( sum( bsxfun(@minus, datawzero ,[X(n,m),Y(n,m)]).^2, 2 ) ); 
    ii = find( dd < len ); 
    cAmat(n,m) = mean(cAdata(ii,:),1); 
end; end 
cAmat = interp2( X', Y', cAmat', XX', YY', 'linear', 1 )'; 


load( 'data_Matlab/DPT_Paul.mat' ) 
dxdynew = dxdy; 
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
