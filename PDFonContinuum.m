% PDF model of HSC differentiation on continuum cell state space 

% Inputs:
%   T: final time 
%   ntest: choose single-cell data set 
%        = 0 - Nestorowa data / 1 - Paul data 
%   nAML : choose AML condition 
%        = INF - normal condition / integer number - AML condition at nAML 
% 
% Outputs:

function PDFonContinuum( T, ntest, nAML ) 
global xx yy N1 N 
%% load scRNA-seq data and Dimension Reduction 

Get_DimReduce_Space( ntest ); 

Compute_systemMatrix( N ); 

uhmsts = Compute_homoestasis( ntest ); 

Compute_systemCoeff( ntest, uhmsts ); 


end 

function Get_DimReduce_Space( ntest ) 
global N 

%%%% default case Nestorowa 
if( nargin == 0 ); ntest = 0; end 

%% load scRNAseq data 
if( ntest == 0 ) 
    filename = 'Nestorowa2016_scRNAseqData.mat'; 
    load(  filename  );    
    DMsig = Inf; 
    N = [100 125]; 
    
elseif( ntest == 1 )
    filename = 'Paul2015_scRNAseqData.mat'; 
    load( filename  );    
    %%%% change data to log2 scale 
    scdata = log2( scdata + 1 );
    DMsig = 50; 
    N = [100 120]; 

end 

%% Dimension Reduction - Diffusion map 
    DMdim = 20; alpha = 0.5; 
   [dc, dcSig] = DiffusionMap_wNewData( scdata', [], DMdim, alpha, DMsig );
    dc = dc_normalize(dc, DMdim);
    
    save( filename, 'dc', 'dcSig' ,'-append' ) 
    
end 

function uhmsts = Compute_homoestasis( ntest ) 
global xx yy N1 

if( ntest == 0 ) 
    filename = 'Nestorowa2016_scRNAseqData.mat';     
elseif( ntest == 1 )
    filename = 'Paul2015_scRNAseqData.mat'; 
end 

 load(  filename  );    
[YY,XX] = meshgrid(yy,xx);
[uhmsts,~,~] = ksdensity( [dc(:,1:2)], [XX(:),YY(:)], 'bandwidth', [0.03 0.03] );
 uhmsts = reshape( uhmsts, N1(1), N1(2) ); 

end 
