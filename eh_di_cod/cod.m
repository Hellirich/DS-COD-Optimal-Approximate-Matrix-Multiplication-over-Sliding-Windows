function [X, Y, idx] = cod(flag, X, Y, l) 
    [mx, n] = size(X);
    my = size(Y, 1);
    
    X = full(X);
    Y = full(Y);   
    
    %%% start QR steps $$$$$$$
    % [Qx, Rx] = qr(X, 0);
    % [Qy, Ry] = qr(Y, 0);
    % we can acclerate QR step if m << d
    [U, S, ~] = svd(X'*X,'econ');
    S(S<1e-10) = 0;
    rk = nnz(S);
    s = diag(S);
    Qx = X * (U(:, 1:rk) * (diag(1./ sqrt(s(1:rk)))));
    Rx = diag(sqrt(s(1:rk))) * U(:, 1:rk)';

    [U, S, ~] = svd(Y'*Y,'econ');
    S(S<1e-10) = 0;
    rk = nnz(S);
    s = diag(S);
    Qy = Y * (U(:, 1:rk) * (diag(1./ sqrt(s(1:rk)))));
    Ry = diag(sqrt(s(1:rk))) * U(:, 1:rk)';
    %%%%%%% end QR steps %%%%%%%%%
    
    [U, S, V] = svd(Rx*Ry','econ');            
    s = diag(S);
    s(s<1e-10) = 0;
    rk = nnz(s);
    if flag == 1
        idx = 2 * l;
        if rk>=l
            X(:, 1:l-1) = Qx * (U(:,1:l-1) * diag(sqrt(s(1:l-1)-s(l))));
            X(:, l:end) = zeros(mx, l+1);                
            Y(:, 1:l-1) = Qy * (V(:,1:l-1) * diag(sqrt(s(1:l-1)-s(l))));
            Y(:, l:end) = zeros(my, l+1);
            idx = l - 1;
        else
            X(:, 1:rk) = Qx * (U(:,1:rk) * diag(sqrt(s(1:rk))));
            X(:, rk+1:end) = zeros(mx, 2*l-rk);                
            Y(:, 1:rk) = Qy * (V(:,1:rk) * diag(sqrt(s(1:rk))));
            Y(:, rk+1:end) = zeros(my, 2*l-rk);   
            idx = rk;
        end
    else
        if rk>=(l+1)
            X(:, 1:l) = Qx * (U(:,1:l) * diag(sqrt(s(1:l)-s(l+1))));
            Y(:, 1:l) = Qy * (V(:,1:l) * diag(sqrt(s(1:l)-s(l+1))));
            X = X(:, 1:l);
            Y = Y(:, 1:l);
        else
            X(:, 1:rk) = Qx * (U(:,1:rk) * diag(sqrt(s(1:rk))));
            Y(:, 1:rk) = Qy * (V(:,1:rk) * diag(sqrt(s(1:rk))));
            X = X(:, 1:rk);
            Y = Y(:, 1:rk);
        end
    end
end
