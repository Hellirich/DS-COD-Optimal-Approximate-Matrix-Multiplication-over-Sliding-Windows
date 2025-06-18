function [A, B, idx, s] = cs(A, B, l)

    [mx,n] = size(A);
    my = size(B, 1);
    [Qx, Rx] = qr(A, 0);
    [Qy, Ry] = qr(B, 0);
    Rx = full(Rx);
    Ry = full(Ry);
    [U, S, V] = svd(Rx*Ry', 'econ');            
    S(S<1e-10) = 0;
    s = diag(S);
    rk = nnz(s);

    if rk>=l
        s = s(1:l)-s(l);
        A(:,1:l) = Qx * (U(:,1:l) * diag(sqrt(s)));
        A(:,l+1:n) = zeros(mx, n-l);
        B(:,1:l) = Qy * (V(:,1:l) * diag(sqrt(s)));
        B(:,l+1:n) = zeros(my, n-l);
        idx = l-1; 
    else
        s = s(1:rk);
        A(:,1:rk) = Qx * U(:,1:rk) * diag(sqrt(s));
        A(:,rk+1:n) = zeros(mx, n-rk);
        B(:,1:rk) = Qy * V(:,1:rk) * diag(sqrt(s));
        B(:,rk+1:n) = zeros(my, n-rk);
        idx = rk;
    end   
end