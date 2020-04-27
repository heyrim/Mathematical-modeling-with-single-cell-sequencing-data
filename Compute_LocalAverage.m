%% Compute local average of cdata on X with distance parameter len 
% Input 
% xdata : two-column vector of 2D location of cdata 
% cdata : vector of function value to interpolate. 
% X     : local average computed on X{1}, X{2} meshgrid 
% len   : distance to compute the local average. 

function cmat = Compute_LocalAverage( xdata, cdata, X, len ) 

cmat = zeros( size(X{1},1), size(X{1},2) ); 
for n = 1:size(X{1},1)
    for m = 1:size(X{1},2) 
        dd = sqrt( sum( bsxfun(@minus, xdata ,[X{1}(n,m),X{2}(n,m)]).^2, 2 ) ); 
        cmat(n,m) = mean(cdata(dd < len,:),1); 
    end
end

end 