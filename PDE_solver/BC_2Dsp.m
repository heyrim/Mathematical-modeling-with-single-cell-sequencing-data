function nC = BC_2Dsp( nC )

nC( [1,end], :) = 0;
if( size( nC,2 ) > 1 )
    nC( :, [1,end]) = 0;
end

end