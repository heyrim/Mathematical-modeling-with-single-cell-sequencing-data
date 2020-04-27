function [cPmat, cAmat, cDth, cVvec, indSig] = PrePocess_Nestorowa(uhmsts)  
global xx yy 

load( 'data/Nestorowa2016_scRNAseqData.mat' ) 
ncluster( ncluster == 9 ) = 2; % ignore small cluster 9 
dc = dc(:,1:2); 


%% Mesh to compute coefficients 
[XX{2},XX{1}] = meshgrid(yy,xx);

%% Compute coefficient on sparser grid (for smoothness) 
len = 0.04; 
[X{2},X{1}] = meshgrid([len/2:(len):(1-len/2)], [len/2:(len):(1-len/2)]); 

%%%% add zeros on coarser grids 
dcwzero = [dc; X{1}(:),X{2}(:)];  
nzero = size(dcwzero,1); 


%% Reaction term - growth rate  
%%%% proliferation rate 
%%%% log(2)/ [8.8*7  12.25  4.3  1] 
cProlif = [0.01125  0.05658  0.1612  0.6931]; 

%%%% Assign rates to each cluster 
cellind2state = [2 1 3 4 4 3 2 3]; 
cPdata = cProlif( cellind2state( ncluster ) )'; 
%%%% append trivial data 
cPdata( (length(cPdata)+1):nzero ) = 0; 

%%%% Local average and interpolation 
cPmat = Compute_LocalAverage( dcwzero, cPdata, X, len ); 
cPmat = interp2( X{1}', X{2}', cPmat', XX{1}', XX{2}', 'linear', 0 )'; 



%% Self renewal fraction 
cAdvec = [0.77  0.7689  0.7359  0.66]; 
cAdata = cAdvec( cellind2state( ncluster ) )'; 
cAdata( (length(cAdata)+1):nzero ) = 1; 

%%%% Local average and interpolation 
cAmat = Compute_LocalAverage( dcwzero, cAdata, X, len ); 
cAmat = interp2( X{1}', X{2}', cAmat', XX{1}', XX{2}', 'linear', 1 )'; 


%% Diffusion pseudotime  - direction of differentiation 
%%%%% reference find_DPF_vector.m 
%%%% load precomputed pseudotime vector 
load('./data/Newtorowa_DiffusionPseudoTime.mat') 

DPTvec(nzero, 2) = 0; 
cVvec = cell(1,2); 
for n = 1:2 
    cVmat(:,:,n) = Compute_LocalAverage( dcwzero, DPTvec(:,n), X, len ); 
    cVvec{n} = interp2( X{1}', X{2}', cVmat(:,:,n)', XX{1}', XX{2}', 'linear', 0 )'; 
    cVvec{n}(isnan(cVvec{n})) = 0;
end 

%% logistic death term 
uhmsts(uhmsts<=eps)=eps;
cDth = 1./uhmsts;

%% cluster index of differentiated cells  
indSig = [3:6,8];
    
end 


