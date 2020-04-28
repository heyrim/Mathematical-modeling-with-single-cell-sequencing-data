function [mappedX, S, mappedY] = DiffusionMap_wNewData(X, Y, no_dims, alpha, sigma)
    
%% Normalize data 
    minX = min(X(:)); 
    X = X - minX; maxX = max(X(:)); 
    X = X / maxX; 
    
if( ~isempty(Y) ) 
    Y = (Y-minX) / maxX; 
end
    
%% Compute kernel matrix 
    
    if( sigma == Inf ) 
    %%%% cosine distance  \eqn{1-corr(c_1, c_2)}   
       [~,m] = size(X); 
        sumX  = sum(X,2)/m; 
        rescaleX = X - sumX*ones(1,m); 
        sumX2 = sqrt(sum( (rescaleX).^2, 2 )/m); 
        K = rescaleX * rescaleX' / m;
        K = K./(sumX2*sumX2');
    
    else    
    %%%% Gaussian distance with sigma 
        sumX = sum(X .^ 2, 2); 
        K = exp(-bsxfun(@plus, sumX, bsxfun(@plus, sumX', -2 * (X * X'))) ./ (2 .* sigma ^ 2)); 

    end 
    
%%%% Similar operation with new data Y 
if( ~isempty(Y) ) 
    if( sigma == Inf ) 
       [~,m] = size(Y); 
        sumY  = sum(Y,2)/m; 
        rescaleY = Y - sumY*ones(1,m); 
        sumY2 = sqrt(sum( (rescaleY).^2, 2 )/m); 
        KK = rescaleY * rescaleX' / m;
        KK = KK./(sumY2*sumX2');
    else
        sumY = sum(Y .^ 2, 2); 
        KK = exp(-bsxfun(@plus, sumY, bsxfun(@plus, sumX', -2 * (Y * X'))) ./ (2 .* sigma ^ 2));     
    end     
end 

%% Compute Markov probability matrix 
    p  = sum(K,1)'; 
    K  = K ./ ((p*p') .^ alpha); 
if( ~isempty(Y) ) 
    q  = sum(KK,2); 
    KK = KK./ ( (q*p').^ alpha ); 
end     
    
%%%%%%  itself probability diag = 0 
%     K = K - diag( diag(K) );  
    p  = sqrt(sum(K,1))'; 
    K  = K ./ (p*p');     
if( ~isempty(Y) ) 
    q  = sqrt(sum(KK,2)); 
    KK = KK./ (q*p'); 
end 

%% Perform SVD 
   [U, S, V] = svd(K, 0); 
    if( size( K, 1 ) > 2500 )
        [U, S, V] = svds(K, no_dims*2); 
    end 

%% Diffusion Components 
    mappedX = bsxfun(@rdivide, U, U(:,1));     
    mappedX = mappedX(:,2:no_dims + 1); 
    
if( ~isempty(Y) ) 
    % operator that gives back DC coordinates : K(1,:) --> U(1,:) 
    % e.g. (K(1,:)*U)./diag(S)' == U(1,:); 
    DCcoeff = (KK*U)./(diag(S))'; % i.e., KK*U = DCcoeff*diag(S); 
    mappedY = bsxfun(@rdivide, DCcoeff, DCcoeff(:,1)); 
    mappedY = mappedY(:,2:no_dims + 1); 
else
    mappedY = []; 
end 

    S = diag(S); S = S(2:no_dims + 1); 
    
    
end 
