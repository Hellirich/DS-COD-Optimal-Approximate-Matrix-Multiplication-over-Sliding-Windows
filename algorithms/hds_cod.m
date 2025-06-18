function result = hds_cod(X, Y, N, l, L_estimate, R_estimate, gap)
    [mx, my, n, ~] = info(X, Y);
    L = ceil(log2(R_estimate));
    LL = max(floor(log2(L_estimate)),1);
    
    err_sum = 0;
    max_err = 0;
    count = 0;
    maxsize = 0;
    updateTime = 0;
    queryTime = 0;
    maxMemoryUsage = 0;

    dscod = cell(1, L);  
    dscod_aux = cell(1, L); 
    theta = cell(1, L);
    for j = 1:L
        theta{j} = 2^(j-1) * (N / l);
    end
    
    for j = LL:L
        ds = FastDumpSnapshotCod(N, mx, my, l, theta{j}, 0);
        ds_aux = FastDumpSnapshotCod(N, mx, my, l, theta{j}, 0);
        dscod{j} = ds;
        dscod_aux{j} = ds_aux;
    end
    
    h = waitbar(0, 'start...'); 
    for i = 1:n
        x = X(:, i);
        y = Y(:, i);
        
        tic;  % update start time
        
        for j = LL:L
            dscod{j}.expire_S(i);
            dscod{j}.cap_S();
            dscod_aux{j}.cap_S();
        end
        
        for j = LL:L
            dscod{j}.update(x, y, i);
            dscod_aux{j}.update(x, y, i);
        end

        if mod(i, N) == 0
            for j = LL:L
                % fprintf("Level %d, S %d, S_aux %d \n", j , dscod{j}.lenS, dscod_aux{j}.lenS)
                dscod{j}.restart(dscod_aux{j})
                dscod_aux{j} = FastDumpSnapshotCod(N, mx, my, l, theta{j}, i);
            end
        end
        
        updateTime = updateTime + toc; % update end time

        % query
        if mod(i, gap) == 0  && i >= N 
           Xw = X(:, max(1, i-N+1):i);
           Yw = Y(:, max(1, i-N+1):i);
           
           tic; % query start time
           for j = LL:L
               if (dscod{j}.lenS < 3*l)
                    break;
               end
           end
           [A,B] = dscod{j}.query();
           queryTime = queryTime + toc; % query end time
           
           XFYF = norm(Xw,'fro') * norm(Yw,'fro');
           err = norm_XYT_ABT(Xw', Yw', A', B')/XFYF;
           % fprintf("S %d, S_aux %d, ",dscod{j}.lenS,dscod_aux{j}.lenS)
           % fprintf('window %d-%d, level %d, l %d, the error is %f\n', max(1, i-N+1), i, j, l, err);
           err_sum = err_sum + err;
           max_err = max(max_err, err);
           count = count + 1;
           sizeW = 0;
           lens = 0;
           for j = LL:L
               lens = lens + dscod{j}.lenS + dscod_aux{j}.lenS;
               sizeW = sizeW + dscod{j}.lenS + dscod_aux{j}.lenS + 2*size(dscod{j}.sketchA, 2) + 2*size(dscod_aux{j}.sketchA, 2) ;
           end
           maxsize = max(maxsize, sizeW);
           
           clear Xw Yw;
           byteStream = getByteStreamFromArray(dscod);
           dscod_mem = length(byteStream);
           byteStream_aux = getByteStreamFromArray(dscod_aux);
           dscod_aux_mem = length(byteStream_aux);
           currentMemory = (dscod_mem + dscod_aux_mem) / 1024 / 1024; 
           if currentMemory > maxMemoryUsage
               maxMemoryUsage = currentMemory;
           end 

        end

        waitbar(i / n, h, sprintf('complete %d/%d (%.2f%%)', i, n, (i / n) * 100));
    end
    avg_err = err_sum / count;
    avgUpdateTime = updateTime / N;
    avgQueryTime = queryTime / count;

    fprintf('hDS-COD, l=%d, avg_error %f, max_error %f, maxsize %d, update %f, query %f, memory %f\n', l, avg_err, max_err, maxsize, avgUpdateTime, avgQueryTime, maxMemoryUsage);
    result = [l, maxsize, avg_err, max_err, maxMemoryUsage, avgUpdateTime, avgQueryTime];
    
    close(h)
end