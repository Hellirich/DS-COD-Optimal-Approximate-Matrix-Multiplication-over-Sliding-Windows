function [L_estimate, R_estimate] = checkR(X, Y)

    norms_X = sqrt(sum(X.^2, 1));   
    norms_Y = sqrt(sum(Y.^2, 1));  
    prodNorm = norms_X .* norms_Y;  
    R_max = max(prodNorm);
    R_avg = mean(prodNorm);
    R_max = full(R_max);
    R_avg = full(R_avg);
    
    fprintf('R_avg = %f R_max = %f\n', R_avg, R_max);
    
    mid = ceil(log2(R_avg));
    high = max(floor(log2(R_max)), mid);

    R_estimate = 2^(floor((mid+high)/2));
    L_estimate = 2^(mid - 2);
    

end