function data = dc_normalize(dc, dim) 

for n = 1:dim 
    minmax(n,1) = min(dc(:,n)); minmax(n,2) = max(dc(:,n)); 
    dL = minmax(n,2) - minmax(n,1);
    minmax(n,1) = minmax(n,1) - dL*0.15; 
    minmax(n,2) = minmax(n,2) + dL*0.15; 
    
end

for n = 1:dim 
    data(:,n) = (dc(:,n)-minmax(n,1))/(minmax(n,2)-minmax(n,1)); 
end


end  