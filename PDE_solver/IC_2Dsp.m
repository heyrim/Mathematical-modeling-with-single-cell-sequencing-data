function uinit = IC_2Dsp( ninit, uHS, ntest )
global xx yy dx N1 N

if( ninit ) 
    
    if( ntest == 0 )
        load( 'data/Nestorowa2016_scRNAseqData.mat' );     
        
       [YY,XX] = meshgrid(yy,xx);
        
        %%%% find Stem cell cluster 
        ii = find( ncluster == 2 );

    elseif( ntest == 1 )
        load( 'data/Paul2015_scRNAseqData.mat' ); 
        
%         y = [-0.2:1.2/N(2):1]';
       [YY,XX] = meshgrid(yy,xx);
        ncluster( ncluster == 10 ) = 6; ii = find( ncluster == 6 );

        
    end
    
    %%%% interpolate Stem cell density
   [uinit,~,~] = ksdensity( [dc(ii,1:2)], [XX(:),YY(:)], 'bandwidth', [0.03 0.03] );
    uinit = reshape( uinit, N1(1), N1(2) );
    
    uinit = ninit * uinit/ (sum(sum(uinit))*prod([dx]));
    uinit = uinit + uHS;
    
else
    
    uinit = uHS;
    
end 


end