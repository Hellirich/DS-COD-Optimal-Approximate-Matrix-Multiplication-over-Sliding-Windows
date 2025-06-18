function ads_cod_time(X, Y, timestamps, N, l, L_estimate, R_estimate, gap)
    lower_F_norm = L_estimate;
    start_time = tic;
    [mx, my, n, ~] = info(X, Y);
    err_sum = 0;
    max_err = 0;
    count = 0;
    maxsize = 0;
    theta = lower_F_norm;

    dscod = FastDumpSnapshotCod(N, mx, my, l, theta, 0);
    dscod_aux = FastDumpSnapshotCod(N, mx, my, l, theta, 0);
    
    t_max = timestamps(end);
    current_data_idx = 1;
    
    for t = 0:t_max

        if mod(t-1, l) == 0
            window_indices = find(timestamps >= t & timestamps <= min(t+l-1,t_max));
            if ~isempty(window_indices)
                fX = X(:, window_indices);
                fY = Y(:, window_indices);
                size_l = norm(fX,'fro') * norm(fY,'fro');
            else
                size_l = 0;
            end
            estimate_L = ceil(log(size_l/l));
            loglogR = log(log(R_estimate));
    
            dscod.nadapt_level(max(estimate_L - loglogR, 0));
            dscod_aux.nadapt_level(max(estimate_L - loglogR, 0));
            
        end

        dscod.expire_S(t);
        
        while current_data_idx <= n && timestamps(current_data_idx) == t
            x = X(:, current_data_idx);
            y = Y(:, current_data_idx);
            dscod.update(x, y, t);
            dscod_aux.update(x, y, t);
            current_data_idx = current_data_idx + 1;
        end

        % N-restart
        if mod(t, N) == 0
            dscod.restart(dscod_aux);
            dscod_aux = FastDumpSnapshotCod(N, mx, my, l, theta, t);
        end
        
        % query
        if mod(t, gap) == 0 && t >= N
            window_start_time = t - N + 1;
            window_indices = find(timestamps >= window_start_time & timestamps <= t);
            if ~isempty(window_indices)
                Xw = X(:, window_indices);
                Yw = Y(:, window_indices);
                [A, B] = dscod.query();
                % err = correlation_error(Xw, Yw, A, B);
                XFYF = norm(Xw, 'fro') * norm(Yw, 'fro');
                err = norm_XYT_ABT(Xw', Yw', A', B')/XFYF;
                % fprintf("S %d, S_aux %d, ",dscod.lenS,dscod_aux.lenS)
                fprintf('window %d-%d, level %d, l %d, the error is %f\n', max(0, t-N+1), t, dscod.level, l, err);
                err_sum = err_sum + err;
                max_err = max(max_err, err);
                count = count + 1;
                sizeW = dscod.lenS + dscod_aux.lenS + 2*size(dscod.sketchA,2) + 2*size(dscod_aux.sketchA,2) ;
                maxsize = max(maxsize, sizeW);
            end
        end
    end
    end_time = toc(start_time);
    avg_err = err_sum / count;
    
    fprintf('aDS-COD_time, for l=%d, the avg_error is %f, the max_error is %f, maxsize is %d, time cost %f\n', ...
        l, avg_err, max_err, maxsize, end_time);
end