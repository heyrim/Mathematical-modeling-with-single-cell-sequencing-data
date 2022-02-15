% PDF model of HSC differentiation on continuum cell state space 

% Inputs:
%   T: final time 
%   ntest: choose single-cell data set 
%        = 0 - Nestorowa data / 1 - Paul data 
%   nAML : choose AML condition 
%        = INF - normal condition / integer number - AML condition at nAML 
% 
% Outputs: 

function PDFonContinuum( N, ntest ) 

close all 
clear all 

global printon cDiff 

%% load scRNA-seq data and Dimension Reduction 

%%%% default case Nestorowa 
addpath('data/', 'dimension_reduction/', 'PDE_solver/', 'misc/'); 

if( nargin < 1 ); N = [100, 125]; end 
if( nargin < 2 ); ntest = 0; end 

if( ntest == 0 ) 
    ll = { 'MPP' 'HSC' 'Baso' 'Ery' 'Neu/Mo'  'Lymph' 'Mk'  'LMPP'  }; 
elseif( ntest == 1 )
    ll = {};     
end


set_system_paramters; 
    
%% Dimension reduction with scRNA data 
Get_DimReduce_Space( ntest ); 

%% Compute PDE system matrix 
Compute_systemMatrix( N ); 

%% Compute homoestasis distribution 
uHS = Compute_homoestasis( ntest ); 

%% Compute PDE coefficient matrix 
Compute_systemCoeff( ntest, uHS ); 

%% IC for normal case - initiate from stem cells 
disp( 'Normal hematopoiesis starting...' ) 

uinit = IC_2Dsp( 0.5, uHS*0.05, ntest ); 

uinit = BC_2Dsp( uinit );

T = 5; dt = 0.001; Tstep = 1;
[uPlot, nCellPlot] = Time_Integ_RK4( @Compute_du_2Dsp, @BC_2Dsp, uinit, T, dt, Tstep, @Compute_cluster ); 

if( printon ) 
figure; hold on;
for nT = 1:(T/Tstep+1) 
    subplot( 1, 5, min(nT, 5));
    imagesc( squeeze(uPlot(:,:,nT))' ); set(gca, 'YDir','normal');  
    title( strcat( 't=', int2str((nT-1)*Tstep) )); 
    xlabel( 'DC1' ); xticks( [] ); ylabel( 'DC2' ); yticks( [] ); 

end
%%%% plot cluster 
% figure; plot( nCellPlot' ); 
% xlabel( 'time' ); ylabel( 'cell number' ); 
% legend( ll ); 
% set( gca, 'fontsize', 14 )

end 


%% AML simulation 
disp( 'AML simulation starting....' )

Update_AML( ntest )

% IC for AML case - initiate from homoestasis
uinit = IC_2Dsp( 0, uHS, ntest ); 

uinit = BC_2Dsp( uinit );

T = 5; dt = 0.001; Tstep = 1;
[uPlot, nCellPlot] = Time_Integ_RK4( @Compute_du_2Dsp, @BC_2Dsp, uinit, T, dt, Tstep, @Compute_cluster ); 

if( printon ) 
figure; hold on;
for nT = 1:(T/Tstep+1) 
    subplot( 1, 5, min(nT, 5));
    title( strcat( 't=', int2str((nT-1)*Tstep) )); 
    imagesc( squeeze(uPlot(:,:,nT))' ); set(gca, 'YDir','normal');  
end
% figure; plot( [0:Tstep:T], nCellPlot' ); 
% xlabel( 'time' ); ylabel( 'cell number' ); 
% legend( ll ); 
% set( gca, 'fontsize', 14 )
end 



%% AML simulation 2 
global cDiff  caml raml 
raml = 0; 
caml = 5; 
cDiff=0.1^3; 
Compute_systemCoeff( ntest, uHS );

Update_AML2( ntest ); 

uinit = IC_2Dsp( 0, uHS, ntest ); 

uinit = BC_2Dsp( uinit );

T = 10; dt = 0.001; Tstep = 2;
[uPlot, nCellPlot] = Time_Integ_RK4( @Compute_du_2Dsp, @BC_2Dsp, uinit, T, dt, Tstep, @Compute_cluster ); 

if( printon ) 
figure; hold on;
for nT = 1:(T/Tstep+1) 
    subplot( 1, 5, min(nT, 5));
    title( strcat( 't=', int2str((nT-1)*Tstep) )); 
    imagesc( squeeze(uPlot(:,:,nT))' ); set(gca, 'YDir','normal');  
end
% figure; plot( [0:Tstep:T], nCellPlot' ); 
% set( gca, 'fontsize', 14 )
% xlabel( 'time' ); ylabel( 'cell number' ); 
% legend( ll ); 

end 
 

%% AML simulation 3 
global cDiff  caml raml 
raml = 1; 
caml = 5; 
cDiff=0.1^3;  
Compute_systemCoeff( ntest, uHS );

Update_AML3( ntest ); 

uinit = IC_2Dsp( 0, uHS, ntest ); 

uinit = BC_2Dsp( uinit );

T = 10; dt = 0.001; Tstep = 2;
[uPlot, nCellPlot] = Time_Integ_RK4( @Compute_du_2Dsp, @BC_2Dsp, uinit, T, dt, Tstep, @Compute_cluster ); 

if( printon ) 
figure; hold on;
for nT = 1:(T/Tstep+1) 
    subplot( 1, 5, min(nT, 5));
    title( strcat( 't=', int2str((nT-1)*Tstep) )); 
    imagesc( squeeze(uPlot(:,:,nT))' ); set(gca, 'YDir','normal');  
end
% figure; plot( [0:Tstep:T], nCellPlot' ); 
% set( gca, 'fontsize', 14 )
% xlabel( 'time' ); ylabel( 'cell number' ); 
% legend( ll); 
end 



end 

function Get_DimReduce_Space( ntest ) 
global N 

%% load scRNAseq data 
    if( ntest == 0 ) 
        filename = 'data/Nestorowa2016_scRNAseqData.mat'; 
        load(  filename  );    
        DMsig = Inf; 
        N = [100 125]; 

    elseif( ntest == 1 )
        filename = 'data/Paul2015_scRNAseqData.mat'; 
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
    filename = 'data/Nestorowa2016_scRNAseqData.mat';     
elseif( ntest == 1 )
    filename = 'data/Paul2015_scRNAseqData.mat'; 
end 

 load(  filename  );    
[YY,XX] = meshgrid(yy,xx);
[uhmsts,~,~] = ksdensity( [dc(:,1:2)], [XX(:),YY(:)], 'bandwidth', [0.03 0.03] );
 uhmsts = reshape( uhmsts, N1(1), N1(2) ); 

end 







