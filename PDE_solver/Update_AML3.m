function Update_AML3( ntest )
global xx yy int_dx cRct cAdv cDth Vaml V V2 addcRct edge N1 cK barv caml raml  

if( ntest == 0 )
    load( 'data/Data_Nestorowa.mat' )
    
    [X,Y] = meshgrid([0:.01:1],[0:.01:1]);
    [XX,YY] = meshgrid(xx,yy);
    MEP = interp2( X,Y, VV{4}, XX, YY )';
        
%     %%%%%%% psuedo inverse
%     load( 'data/Newtorowa_AML_vec.mat' );
    
    %%%%%%% potential
    % % global cDiff
    AMLcenter = [0.6101   0.2153]; % shiifted from [1 -2]    
    usol = Normal2D(AMLcenter, 0.01*[1,0.3;0.3,1], xx', yy' );
    
    Vaml = Compute_V( usol ) ;
    
%     V{1}(:) = 0; V{2}(:) = 0; 
    V2{1}(:) = 0; V2{2}(:) = 0; 
	
    usol(usol>4) = 4;
    addcRct = raml * usol/max(usol(:)) * 0.6931;  
    cDth = cDth .* (1-usol/max(usol(:))); 
    
elseif( ntest == 1 )
    
    load( 'data/Data_Paul.mat' );
    MEP(:, 21:121) = WW;
    %  STHSC(:, 21:121) = VV{6};
    
    %  load( 'data_Matlab/190604_Paul_AML_vec.mat' );
    %  load( 'data_Matlab/190612_Paul_AML_vec_straight.mat' );
    load( 'data/190614_Paul_AML_vec.mat' );
    for n = 1:2; Vaml{n}(:,21:121) = newvec(:,:,n); end 
    
    %%%%%% potential
    AMLcenter = [0.6  0.2];
    usol = Normal2D(AMLcenter, 0.04*[1,0;0,0.2], xx', yy' );
    %%%%%% vector from potential
    Vaml = Compute_V( usol ) ;
    
    
    usol(usol>4) = 4;
    addcRct = raml * usol/max(usol(:)) * 0.6931;
    cDth = cDth .* (1-usol/max(usol(:)));     
    
end


for n = 1:2; Vaml{n} = caml *Vaml{n}; end



end