function uinit = IC_2Dsp( ninit, uHS, ntest )
global xx yy dx N1 N

if( ninit ) 
    
    if( ntest == 0 )
        load( 'data/Nestorowa2016_scRNAseqData.mat' );     

        [X,Y] = meshgrid([0:.01:1],[0:.01:1]);
        [YY,XX] = meshgrid(yy,xx);
        ii = find( ncluster == 2 );

    elseif( ntest == 1 )
        load( 'data/Paul2015_scRNAseqData.mat' ); %dc(:,2) = (dc(:,2)+0.2)*(5/6);
        
        y = [-0.2:1.2/N(2):1]';
        [YY,XX] = meshgrid(y,xx);
        ncluster( ncluster == 10 ) = 6; ii = find( ncluster == 6 );

        
    end
    

   [uinit,~,~] = ksdensity( [dc(ii,1:2)], [XX(:),YY(:)], 'bandwidth', [0.03 0.03] );
    uinit = reshape( uinit, N1(1), N1(2) );
    
    uinit = ninit * uinit/ (sum(sum(uinit))*prod([dx]));
    uinit = uinit + uHS;
    
else
    
    uinit = uHS;
    
end 


end