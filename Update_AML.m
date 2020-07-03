function indSig = Update_AML( ntest )
global xx yy int_dx cRct cAdv cDth Vaml V addcRct edge N1 cK barv caml

if( ntest == 0 )
    load( 'data_Matlab/181019_Data_Nestorowa.mat' )
    
    [X,Y] = meshgrid([0:.01:1],[0:.01:1]);
    [XX,YY] = meshgrid(xx,yy);
    MEP = interp2( X,Y, VV{4}, XX, YY )';
    STHSC = interp2( X,Y, VV{1}, XX, YY )';
    
    %%%%%%% signal param
    % barv = barv*10;
    % cK = 10;
    
    %%%%%%% psuedo inverse
    load( 'data_Matlab/190406_Nstrw_AML_vec.mat' );
    
    for n = 1:2; Vaml{n} = newvec(:,:,n); end
    
    %%%%%%% potential
%     % % global cDiff
%     AMLcenter = [0.6101   0.2153]; % shiifted from [1 -2]    
%     usol = Normal2D(AMLcenter, 0.04*[1,0.3;0.3,1], xx', yy' );
%     
%     Vaml = Compute_V( usol, 0 ) ;
%     
%     usol(usol>10) = 10;
%     addcRct = usol/max(usol(:)) * 0.6931;
    
    
%      STHSC(STHSC>20)=20; STHSC = STHSC/20 * (1-1/2.5);
%      cDth = cDth.*(1 - STHSC);
    STHSC(MEP>20)=20; STHSC = STHSC/20 * 4 + 1;
    cDth = cDth./STHSC;
    
    
elseif( ntest == 1 )
    
    load( 'data_Matlab/190604_Data_Paul.mat' );
    MEP(:, 21:121) = WW;
    STHSC = zeros( size(WW) ); 
    
    %  load( 'data_Matlab/190604_Paul_AML_vec.mat' );
    %  load( 'data_Matlab/190612_Paul_AML_vec_straight.mat' );
    load( 'data_Matlab/190614_Paul_AML_vec.mat' );
    for n = 1:2; Vaml{n}(:,21:121) = newvec(:,:,n); end
    
    %%%%%%% potential
%     AMLcenter = [0.6  0.2];
%     usol = Normal2D(AMLcenter, 0.04*[1,0;0,0.2], xx', yy' );
%     %%%%%% vector from potential
%     Vaml = Compute_V( usol, 0 ) ;
%     
%     
%     usol(usol>10) = 10;
%     addcRct = usol/max(usol(:)) * 0.6931;
    
    
end

%  MEP(MEP>20)=20; MEP = MEP/20 * (1-0.1);
%  cDth = cDth.*(1 - MEP);
MEP(MEP>20)=20; MEP = MEP/20 * 19 + 1;
cDth = cDth./MEP;



for n = 1:2; Vaml{n} = caml *Vaml{n}; end

for n = 1:2; Vaml{n} = Vaml{n} - V{n}; end
% % ind = find( usol >= 10 );
% % for n = 1:2; Vaml{n}(ind) = Vaml{n}(ind) - V{n}(ind); end


if( ntest == 0 )
    indSig = [3,5,6,8];
elseif( ntest == 1 )
    indSig = [2,3,9];
end

end