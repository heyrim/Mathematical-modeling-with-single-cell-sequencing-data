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

uinit = IC_2Dsp( 0.1, uHS*0.05, ntest ); 
uinit = BC_2Dsp( uinit );

Time_Integ_RK4( @Compute_du_2Dsp, @BC_2Dsp, uinit, T, dt, Tstep ) 




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

function dnC = Compute_du_2Dsp( nC, time )
global Dxn Dxxn cDiff V Dxm Dxxm cRct cDth cAdv cA Vaml cS addcRct

React = cRct .* nC ;
Advec = cS * 2*( 1-cAdv).* React;

React = ( 1 - cDth.*nC ); React(React<-0.6925)=-0.6925;
React = cRct .*React.*nC;

dnC = cDiff *(Dxxn *nC + nC* Dxxm') ...
        - Dxn * ( V{1}.* nC + cA{1}.*Advec) - ( V{2}.* nC + cA{2}.*Advec)*Dxm' + React;

React = addcRct.*nC;
dnC = dnC - Dxn * ( Vaml{1}.*React) - ( Vaml{2}.*React) * Dxm' ; % + React;
% dnC = dnC + React;


end

function nC = BC_2Dsp( nC )

nC( [1,end], :) = 0;
if( size( nC,2 ) > 1 )
    nC( :, [1,end]) = 0;
end

end



function uinit = IC_2Dsp( ninit, uHS, ntest )
global xx yy dxn dxm N1 N

if( ntest == 0 )
    load( 'Nestorowa2016_scRNAseqData.mat' );     

    [X,Y] = meshgrid([0:.01:1],[0:.01:1]);
    [XX,YY] = meshgrid(xx,yy);
    uinit = interp2( X,Y, VV{2}, XX, YY )';
            
elseif( ntest == 1 )
    load( 'Paul2015_scRNAseqData.mat' ); %dc(:,2) = (dc(:,2)+0.2)*(5/6);
    y = [-0.2:1.2/N(2):1]';
    [YY,XX] = meshgrid(y,xx);
    cluster(cluster==10) = 6; ii = find( cluster == 6 );
    [uinit,~,~] = ksdensity( [dc(ii,:)], [XX(:),YY(:)], 'bandwidth', [0.03 0.03] );
    
    %     ii = setdiff( ii, find( dc(:,1) > 0.2 ) ); ii = setdiff( ii, find( dc(:,2) > 0.2 ) );
    %    [uinit,~,~] = ksdensity( [dc(ii,:)], [XX(:),YY(:)], 'bandwidth', [0.02 0.02] );
    
    uinit = reshape( uinit, N1(1), N1(2) );
    
    %     load( 'data_Matlab/190517_PaulInit.mat' );
end

%     nCell = uSH*ninit;
if( nargin == 1 )
    ninit = 1;
end

nCell = uinit/ (sum(sum(uinit))*prod([dxn dxm])) *ninit ;

if( ~isempty(uSH) )
    nCell = nCell + uSH;
end


end

