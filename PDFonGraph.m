%% PDF on graph with multiple end nodes 
function PDFonGraph( Tin, Nin, nfig, ntestin )

close all 
clear all 

global dx BC printon cDiff cAdvM cK cS barv 
    
    addpath( './PDE_solver/', './data/', './dimension_reduction/', './misc/') 

    set_system_paramters_graph;    
    
    if( nargin < 1 ) ;        T = 2;       else; T = Tin; end 
    if( nargin < 2 ) ;        N = 100;     else; N = Nin; end 
    if( nargin < 3 ) ;        nfig = floor( rand(1)*1000 ); end         
    
    dt = 0.1^4;     
    
%% Compute PDE system matrix 
    Compute_systemMatrix_Graph( N );     
    
%% IC & parameters 
    if( nargin > 3 ); ntest = ntestin; else;    ntest = 0;  end %%% ntest == 0 for Nestrowa graph 9 / ntest == 1 for Paul graph 
   [BC, uHS, indSig] = edge_graph( ntest );  

   [uinit] = IC_graph( 0.5, uHS*0.05, ntest );  
        
   [nTotVrtx(:,1), nTotCell(1)] = Compute_nCell_graph( uinit ); 
    
    
    cS = 1-1/(1+exp( -(sum(nTotVrtx(indSig, 1)) - barv)*cK )); 
    cSsave(1) = cS;     

%% Time integration 
    if( dx/max(max(cAdvM))/2 < dt || dx^2/max(cDiff)/2 < dt ) 
        disp( 'reduce dt' ); 
    end 
    
    disp( strcat( 'params : ', num2str([cDiff cS cK barv]) ) ); 
        
    T = 5; dt = 0.0001; Tstep = 1;
    [uPlot] = Time_Integ_RK4( @Compute_du_graph, @BC_graph, uinit, T, dt, Tstep, @Compute_nCell_graph );     

    if( printon ) 
      for nn = 1:size( uPlot, 3 ) 
        plot_graph3d( uPlot(:,:,nn), nfig  , [1,size( uPlot, 3 ),nn] ); caxis( [0 1] )
        title( strcat( 't=', int2str( (nn-1)*Tstep ) )); 
        xlabel( 'DC1' ); xticks( [] ); ylabel( 'DC2' ); yticks( [] ); 
        plot_graph2d( uPlot(:,:,nn), nfig+1, [1,size( uPlot, 3 ),nn] ); caxis( [0 1] )
        title( strcat( 't=', int2str( (nn-1)*Tstep ) )); 
        xlabel( 'x' ); 
      end
    end
    
global nAML 
    if( nAML ) 
        [indSig] = Update_AML( U, ntest ); 
    end 
    

end 

function Update_adv( U, indSig )
global cK cS barv wUn int_dx % cSelfRnw cRct cAdv cAdvM N1 nEdge  

    Un = Compute_nCell_graph( U );  
    cS = 1-1/(1+exp( -(sum(Un(indSig)) - barv)*cK )); 

end 

function indSig = Update_AML( U, ntest )
global int_dx cRctM cAdvM cDthM addcRct cRct cAdv edge vertex N1 cK barv cDiffM cDifff


if( ntest == 0 ) %%%% Graph 8
    indSig = [3,5,6,8];
    
    %%%% Ery 4 / Mk 8 flow block
    ii = find( edge(:,2) == 4 )';
    for n = ii
        cAdvM(:,n) = 0;
        
        cDthM(:,n) = cDthM(:,n)./linspace(1, 20, N1)';
        
    end
    %%%%% ST-HSC
    ii = find( edge(:,2) == 1 )';
    for n = ii
        
        cDthM(:,n) = cDthM(:,n)./linspace(1, 5, N1)';
    end
    
elseif( ntest == 21 ) %%%% Graph 21
    indSig = [19, 9, 6, 10, 18, 8, 17, 13, 5 ];
    
    ii = [find( edge(:,2) == 16 )';find( edge(:,2) == 14 )'];
    for n = ii
        cAdvM(:,n) = 0;
        cDthM(:,n) = cDthM(:,n)./linspace(15, 20, N1)';
    end
    
    ii = find( edge(:,2) == 7 )';
    for n = ii
        cAdvM(:,n) = 0;
        cDthM(:,n) = cDthM(:,n)./linspace(1, 15, N1)';
    end
    
    ii = find( edge(:,2) == 3 )'; %%%%% ST-HSC
    for n = ii
        cDthM(:,n) = cDthM(:,n)./linspace(1, 5, N1)';
    end
    ii = find( edge(:,1) == 3 )'; %%%%% ST-HSC
    for n = ii
        cDthM(:,n) = cDthM(:,n)./linspace(5, 1, N1)';
    end
    
elseif( ntest == 109 )
    
    indSig = [2,3,9];
    %%%% Ery 8, 4 flow block
    ii = setdiff( [find( edge(:,2) == 4 )',find( edge(:,2) == 8 )'], 9 );
    for n = ii
        cAdvM(:,n) = 0;
        cDthM(:,n) = cDthM(:,n)./linspace(1, 20, N1)';
    end
    cAdvM(:,9) = 0;
    cDthM(:,9) = cDthM(:,9)./20; %%%%% edge = [4,8]
    
end


% reduce death rate around pro-Meg/E 
%% similar to graph simulation 
%     load( 'data/Nestorowa_G8.mat' ); 
%     vertex = centerdm(1:8,:); 
%     
%     ii = find( edge(:,2) == 4 )'; 
% for n = ii 
% % %     cRctM(:,n) = cRctM(:,n) + tmp; 
% 
% %     tmp = linspace( 0, -0.5, ((N1+1)/2) ); 
% %     tmp = [tmp(1:end-1), linspace( -0.5, 0, ((N1+1)/2) ) ];
% %     cDthM(:,n) = cDthM(:,n) + tmp'; 
%     
% %     tmp = interp1([0,0.5,1], [1 2 1], [0:0.02:1], 'spline' ); 
%     tmp = interp1([0,0.5,1], [0 0.5 0], [0:0.02:1], 'spline' ); 
%     cDthM(:,n) = cDthM(:,n) - tmp'; 
% end 
% ii = find( edge(:,2) == 5 )'; 
% for n = ii 
%     tmp = linspace(  0.25,  0.1, (N1) ); 
%     cRctM(:,n) = cRctM(:,n) + tmp'; 
%     cDthM(:,n) = cDthM(:,n) + tmp'; 
% end 

% barv = barv*10; 
% cK = 10; 

end 
