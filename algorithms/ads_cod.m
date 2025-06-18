function result = ads_cod(X, Y, N, l, R_estimate, gap)
    [mx, my, n, ~] = info(X, Y);
    err_sum = 0;
    max_err = 0;
    count = 0;
    maxsize = 0;
    theta = N/l;
    updateTime = 0;
    queryTime = 0;
    maxMemoryUsage = 0;
    dscod = FastDumpSnapshotCod(N, mx, my, l, theta, 0);
    dscod_aux = FastDumpSnapshotCod(N, mx, my, l, theta, 0);
    
    h = waitbar(0, 'start...'); 
    for i = 1:n
        
        tic;  % update start time
        if mod(i-1, l) == 0
            fX = X(:, i:min(i+l-1,n));
            fY = Y(:, i:min(i+l-1,n));
            size_l = norm(fX,'fro') * norm(fY,'fro');
            estimate_L = ceil(log(size_l/l));
            loglogR = log(log(R_estimate));
    
            dscod.adapt_level(max(estimate_L - loglogR, 0));
            dscod_aux.adapt_level(max(estimate_L - loglogR, 0));
        end
        
        x = X(:, i);
        y = Y(:, i);
        
        dscod.expire_S(i);

        dscod.update(x, y, i);
        dscod_aux.update(x, y, i);
        
        % N-restart
        if mod(i, N) == 0
            dscod.restart(dscod_aux)
            dscod_aux = FastDumpSnapshotCod(N, mx, my, l, theta, i);
        end

        updateTime = updateTime + toc; % update end time

        % query
        if mod(i, gap) == 0 && i >= N
           Xw = X(:, max(1, i-N+1):i);
           Yw = Y(:, max(1, i-N+1):i);

           tic; % query start time
           [A,B] = dscod.query();
           
           queryTime = queryTime + toc; % query end time

           XFYF = norm(Xw, 'fro') * norm(Yw, 'fro');
           err = norm_XYT_ABT(Xw', Yw', A', B')/XFYF;
           % fprintf("S %d, S_aux %d, ",dscod.lenS,dscod_aux.lenS)
           % fprintf('window %d-%d, level %d, l %d, the error is %f\n', max(1, i-N+1), i, dscod.level, l, err);
           err_sum = err_sum + err;
           max_err = max(max_err, err);
           count = count + 1;
           sizeW = dscod.lenS + dscod_aux.lenS + 2*size(dscod.sketchA,2) + 2*size(dscod_aux.sketchA,2) ;
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
    
    fprintf('aDS-COD, l=%d, avg_error %f, max_error %f, maxsize %d, update %f, query %f, memory %f\n', l, avg_err, max_err, maxsize, avgUpdateTime, avgQueryTime, maxMemoryUsage);
    result = [l, maxsize, avg_err, max_err, maxMemoryUsage, avgUpdateTime, avgQueryTime];
    close(h)
end