function value = Normal( mu, sig, x )
    value = 1.0 / sqrt( 2.0*pi ) /sig * exp( -(( x - mu ).^2) /2.0 /sig^2 ) ;
end 