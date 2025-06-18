function [mx, my, n, R_max] = info(X, Y)
    [mx, n] = size(X);
    my = size(Y, 1);
    norms_X = sqrt(sum(X.^2, 1));   
    norms_Y = sqrt(sum(Y.^2, 1)); 
    prodNorm = norms_X .* norms_Y;  
    R_max = max(prodNorm);
    R_max = full(R_max);
end