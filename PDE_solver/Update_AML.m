function Update_AML( ntest )
global xx yy indSig cDth V V2 

if( ntest == 0 )
    load( 'data/Data_Nestorowa.mat' )
    
    [X,Y] = meshgrid([0:.01:1],[0:.01:1]);
    [XX,YY] = meshgrid(xx,yy);
    
    %%%% cluster density 
    MEP = interp2( X,Y, VV{4}, XX, YY )';
    STHSC = interp2( X,Y, VV{1}, XX, YY )';
    
    %%%% exclude MEP from homeostasis 
    indSig = setdiff( indSig, 4 );    % [3,5,6,8];
    
    
    
elseif( ntest == 1 )
    
    load( 'data/Data_Paul.mat' );
    MEP = WW; %% MEP cluster 4 and 8 together 
    
    %%%% exclude MEP from homeostasis 
    indSig = setdiff( indSig, [4,8] );    % [2,3,9];   
    
end

%% Alter growth term of MEP for increased population 
MEP(MEP>20)=20; MEP = MEP/20 * 19 + 1;
cDth = cDth./MEP;

if( ntest == 0 )
    %%%% Alter growth term of STHSC for increased population 
    STHSC(STHSC>20)=20; STHSC = STHSC/20 * 4 + 1;
    cDth = cDth./STHSC;
end 
    
%% Alter advection term for blocked differentiation to MEP 
for n = 1:2 
    V{n}( MEP>20-eps ) = 0; 
    V2{n}( MEP>20-eps ) = 0; 
end 

end