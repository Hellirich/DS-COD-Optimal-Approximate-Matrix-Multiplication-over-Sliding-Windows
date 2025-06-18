function hds_cod_time(X, Y, timestamps, N, l, R_lower_estimated, R_upper_estimated, gap)
    lower_F_norm = R_lower_estimated;    
    upper_F_norm = R_upper_estimated;
    
    start_time = tic;
    err_sum = 0;
    max_err = 0;
    count = 0;
    maxsize = 0;
    [mx, my, n, R] = info(X, Y);
    L = floor(log2(N*R/l));
    base = 1;
    if upper_F_norm ~= 0
        L = floor(log2(upper_F_norm/lower_F_norm));
        base = N*lower_F_norm/l;
    end
    
    
    dscod = cell(1, L+1);
    dscod_aux = cell(1, L+1);
    theta = cell(1, L+1);
    for j = 1:L+1
        theta{j} = base*2^(j-1);
        dscod{j} = FastDumpSnapshotCod(N, mx, my, l, theta{j}, 0);
        dscod_aux{j} = FastDumpSnapshotCod(N, mx, my, l, theta{j}, 0);
    end

    t_max = timestamps(end);
    current_data_idx = 1;
    
    for t = 0:t_max
        for j = 1:L+1
            dscod{j}.expire_S(t);
            dscod{j}.cap_S();
            dscod_aux{j}.cap_S();
        end

        while current_data_idx <= n && timestamps(current_data_idx) == t
            x = X(:, current_data_idx);
            y = Y(:, current_data_idx);
            for j = 1:L+1
                dscod{j}.update(x, y, t);
                dscod_aux{j}.update(x, y, t);
            end
            current_data_idx = current_data_idx + 1;
        end
        
        % N-restart
        if mod(t, N) == 0
            for j = 1:L+1
                % fprintf("Level %d, S %d, S_aux %d \n", j , dscod{j}.lenS, dscod_aux{j}.lenS)
                dscod{j}.restart(dscod_aux{j})
                dscod_aux{j} = FastDumpSnapshotCod(N, mx, my, l, theta{j}, t);
            end
        end

        % query
        if mod(t, gap) == 0  && t >= N 
           window_start_time = t - N + 1;
           window_indices = find(timestamps >= window_start_time & timestamps <= t);
           if ~isempty(window_indices)
               Xw = X(:, window_indices);
               Yw = Y(:, window_indices);
               for j = L+1:-1:1
                   if (dscod{j}.lenS >= 2*l) || (dscod{j}.lenS() >= l && dscod{max(j-1,1)}.lenS() >= 3 * l)
                        break;
                   end
               end
               [A,B] = dscod{j}.query();
               % err = correlation_error(Xw, Yw, A, B);
               XFYF = norm(Xw, 'fro') * norm(Yw, 'fro');
               err = norm_XYT_ABT(Xw', Yw', A', B')/XFYF;
               % fprintf("S %d, S_aux %d, ",dscod{j}.lenS,dscod_aux{j}.lenS)
               % fprintf('window %d-%d, level %d, l %d, the error is %f\n', max(1, t-N+1), t, j, l, err);
               err_sum = err_sum + err;
               max_err = max(max_err, err);
               count = count + 1;
               sizeW = 0;
               for j = 1:L+1
                    sizeW = sizeW + dscod{j}.lenS + dscod_aux{j}.lenS + 2*size(dscod{j}.sketchA, 2) + 2*size(dscod_aux{j}.sketchA, 2) ;
               end
               maxsize = max(maxsize, sizeW);
           end
        end
    end
    avg_err = err_sum / count;
    end_time = toc(start_time);
    fprintf('hDS-COD, for l=%d, the avg_error is %f, the max_error is %f, maxsize is %d, time cost %f\n', l, avg_err, max_err, maxsize, end_time);
end