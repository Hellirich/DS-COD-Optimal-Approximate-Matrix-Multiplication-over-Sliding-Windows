function sigma_1 = norm_XYT_ABT(X, Y, A, B, T)
    if nargin < 5
        T = 100;
    end
        
    v = randn(size(Y, 2), 1);
    v = v / norm(v);
    
    for i=1: T
        v = X'*(Y*v) - A'*(B*v); 
        v = Y'*(X*v) - B'*(A*v); 
        v = v / norm(v);
    end
    
    sigma_1 = norm(X'*(Y*v) - A'*(B*v));
end