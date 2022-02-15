function Update_AML2( ntest )
global xx yy int_dx cRct cAdv cDth Vaml V V2 addcRct edge N1 cK barv caml raml 

if( ntest == 0 )
    load( 'data/Data_Nestorowa.mat' )
    
    [X,Y] = meshgrid([0:.01:1],[0:.01:1]);
    [XX,YY] = meshgrid(xx,yy);
    MEP = interp2( X,Y, VV{4}, XX, YY )';
        
    load( 'data/Newtorowa_AML_vec.mat' );
    
elseif( ntest == 1 )
    
    load( 'data/Data_Paul.mat' );
    [XX,YY] = meshgrid(xx,yy);
    MEP = WW;
    
    load( 'data/Paul_AML_vec.mat' );    
    
end

% V{1}(:) = 0; V{2}(:) = 0; 
V2{1}(:) = 0; V2{2}(:) = 0; 
for n = 1:2; Vaml{n} = newvec(:,:,n); end   
for n = 1:2; Vaml{n} = caml *Vaml{n}; end  


% % allow cells to grow near new shifted location 
ind = find( MEP >= 20-eps ); 
XXind = XX(ind) + Vaml{1}(ind); 
YYind = YY(ind) + Vaml{2}(ind); 

[usol,~,~] = ksdensity( [XXind,YYind], [XX(:),YY(:)], 'bandwidth', [0.02 0.02] );
usol = reshape( usol, N1(1), N1(2) );
usol(usol>20) = 20; usol = usol/20; 
cDth = cDth .* (1-usol); 
addcRct = raml *usol*0.6931; 


end