function result = di_cod(X, Y, N, L_number, gap)
    [mx, my, n, R] = info(X, Y);
    l_avg = 0;
    err_avg = 0;
    err_max = 0;
    l_sum = 0;
    maxsize = -1;
    cnt = 0;
    L = {};
    L_activate = {};
    updateTime = 0;
    queryTime = 0;
    maxMemoryUsage = 0;

    % Initialise each layer
    for i = 1:L_number 
        L_activate{end + 1} = Block(mx, my, 2^(i+1), 0);
        L{end + 1} = {};
    end
    X_size = 0;
    Y_size = 0;
    len = 0;
    
    h = waitbar(0, 'start...'); 
    for t = 1:n
        x = X(:, t);
        y = Y(:, t);

        tic;  % update start time
        for i = 1:L_number
            l = 2^i;
            if ~isempty(L{i}) && L{i}{1}.t_start <= t - N
                l_sum = l_sum - l;
                L{i}(1) =  [];
            end
            L_activate{i}.idx = L_activate{i}.idx + 1;
            L_activate{i}.X(:, (L_activate{i}.idx)) = x;
            L_activate{i}.Y(:, (L_activate{i}.idx)) = y;
            if L_activate{i}.idx >= 2 * l || t == n
                [L_activate{i}.X, L_activate{i}.Y, L_activate{i}.idx] = cod(1, L_activate{i}.X, L_activate{i}.Y, l);
            end 
            L_activate{i}.t_end = t;
        end
        X_size = X_size + norm(x, 2) ^ 2;
        Y_size = Y_size + norm(y, 2) ^ 2;
        
        if X_size > N * R / 2^L_number || Y_size > N * R / 2^L_number || t == n
            X_size = 0;
            Y_size = 0;
            len = len + 1;
            v = min(tail_zeros(len) + 1, L_number);
            for i = 1:v
                l = 2^i;
                [L_activate{i}.X, L_activate{i}.Y, L_activate{i}.idx] = cod(1, L_activate{i}.X, L_activate{i}.Y, l);
                block = Block(mx, my, l, t); 
                block.X = L_activate{i}.X(:, 1:l);
                block.Y = L_activate{i}.Y(:, 1:l);
                block.t_start = L_activate{i}.t_start;
                l_sum = l_sum + l;
                L{i}{end + 1} = block;
                L_activate{i}.clear(t); 
                L_activate{i}.t_start = t;
                L_activate{i}.t_end = t;
                L_activate{i}.X = zeros(L_activate{i}.mx, L_activate{i}.l);
                L_activate{i}.Y = zeros(L_activate{i}.my, L_activate{i}.l);
                L_activate{i}.idx = 0;
            end
        end
        updateTime = updateTime + toc; % update end time

        if t >= N && mod(t, gap) == 0
            Xw = X(:, max(1, t - N + 1):t);
            Yw = Y(:, max(1, t - N + 1):t);
           
            tic; % query start time
            A = [];
            B = [];
            a = t;
            b = t - N;
            v = length(L);
            
            while v > 0 && isempty(L{v})
                v = v - 1;
            end
            
            st = 100000;
            en = -1;
            for i = v:-1:1
                if L{i}{1}.t_end <= a
                    A = horzcat(L{i}{1}.X, A);
                    B = horzcat(L{i}{1}.Y, B);
                    a = L{i}{1}.t_start;
                    st = min(st, a);
                end
                if length(L{i}) == 1
                    b = L{i}{1}.t_end;
                    en = max(en, b);
                elseif L{i}{end}.t_start >= b
                    A = horzcat(A, L{i}{end}.X);
                    B = horzcat(B, L{i}{end}.Y);
                    b = L{i}{end}.t_end;
                    en = max(en, b);
                end
            end
            queryTime = queryTime + toc; % query end time

            XFYF = norm(Xw, 'fro') * norm(Yw, 'fro');
            err = norm_XYT_ABT(Xw', Yw', A', B')/XFYF;
            [~, ll] = size(A);
            sizeW = 0;
            for i = 1:length(L)
                currentLayer = L{i};
                for j = 1:length(currentLayer)
                    element = currentLayer{j};
                    sizeW = sizeW + 2 * size(element.X, 2);
                end
                sizeW = sizeW + 2 * size(L_activate{i}.X, 2); % active buffer size
            end
            
            l_avg = l_avg + ll;
            err_avg = err_avg + err;
            err_max = max(err_max,err);
            cnt = cnt + 1;
            maxsize = max(maxsize, sizeW);

            clear Xw Yw;
            
            byteStream = getByteStreamFromArray(L);
            L_mem = length(byteStream);
            byteStream = getByteStreamFromArray(L_activate);
            L_active_mem = length(byteStream);

            currentMemory = (L_mem + L_active_mem) / 1024 / 1024;
            if currentMemory > maxMemoryUsage
                maxMemoryUsage = currentMemory;
            end
            
        end

        waitbar(t / n, h, sprintf('complete %d/%d (%.2f%%)', t, n, (t / n) * 100));


    end

    l_avg = floor(l_avg / cnt);
    err_avg = err_avg / cnt;
    avgUpdateTime = updateTime / N;
    avgQueryTime = queryTime / cnt;


    fprintf('DI-COD, l=%d, avg_error %f, max_error %f, maxsize %d, update %f, query %f, memory %f\n', l_avg, err_avg, err_max, maxsize, avgUpdateTime, avgQueryTime, maxMemoryUsage);
    result = [l_avg, maxsize, err_avg, err_max, maxMemoryUsage, avgUpdateTime, avgQueryTime];
    close(h)
    
    
end

function idx = tail_zeros(n)
    idx = 0;
    while bitget(n, 1) == 0
        idx = idx + 1;
        n = bitshift(n, -1);
    end
end
