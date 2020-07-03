
function data = dc_normalize(dc, dim) 

    minmax = zeros( dim, 2 ); 
    for n = 1:dim 
        minmax(n,1) = min(dc(:,n)); minmax(n,2) = max(dc(:,n)); 
        dL = minmax(n,2) - minmax(n,1);
        minmax(n,1) = minmax(n,1) - dL*0.15; 
        minmax(n,2) = minmax(n,2) + dL*0.15; 

    end

    data = zeros( size( dc ) ); 
    for n = 1:dim 
        data(:,n) = (dc(:,n)-minmax(n,1))/(minmax(n,2)-minmax(n,1)); 
    end


end  
