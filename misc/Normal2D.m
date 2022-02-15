function value = Normal2D( mu, covC, x, y ) 
    
    detC = det(covC);     invC = inv(covC); 
    
    [xx,yy] = ndgrid( x, y ); xx = xx - mu(1); yy = yy - mu(2); 
    value = 1.0 / ( 2.0*pi ) / sqrt(detC)  ... 
               * exp( -(xx .* invC(1,1) .* xx) /2.0 ) .* exp( -(xx .* invC(1,2) .* yy) /2.0 )  ...
                .* exp( -(yy .* invC(2,1) .* xx) /2.0 ).* exp( -(yy .* invC(2,2) .* yy) /2.0 ) ; 

end 
