function result = eh_cod(X, Y, N, l, gap)
    [mx, n] = size(X);
    my = size(Y, 1);
    maxsize = -1;
    err_avg = 0;
    err_max =0;
    b = l/2;
    q = 0;
    L = {};
    lv1 = Level();
    L{end+1} = lv1;
    active = Block(mx, my, l, 0);
    updateTime = 0;
    queryTime = 0;
    maxMemoryUsage = 0;
    
    h = waitbar(0, 'start...'); 
    for t = 1:n
        x = X(:, t);
        y = Y(:, t);
    
        tic;  % update start time
        % Pop expired blocks
        if  ~isempty(L{end}.blocks) && L{end}.blocks{1}.t_end < t - N
            L{end}.blocks(1) = [];
            if isempty(L{end}.blocks)
                L(end) = [];
            end
        end

        active = active.insert(x, y, t);
    
        if active.size >= l || t == n
            block = Block(mx, my, l, t);
            block.X = active.X(:, 1:active.idx);
            block.Y = active.Y(:, 1:active.idx);
            block.t_start = active.t_start;
            block.sizex = active.sizex;
            block.sizey = active.sizey;
            block.size = active.size;
            L{1}.blocks{end+1} = block;
            active = active.clear(t);
            for i = 1:length(L)
                if length(L{i}.blocks) > b
                    if i + 1 == length(L) + 1
                        new_level = Level();
                        L{end+1} = new_level;
                    end
                    if L{i}.blocks{1}.size > (2^(i))*l
                        block0 = L{i}.blocks{1};
                        L{i}.blocks(1) = [];
                        L{i + 1}.blocks{end+1} = block0;
                    else
                        new_block = Block(mx, my, l, t);
                        new_block.t_start = L{i}.blocks{1}.t_start;
                        new_block.t_end = L{i}.blocks{2}.t_end;
                        new_block.sizex = L{i}.blocks{1}.sizex + L{i}.blocks{2}.sizex;
                        new_block.sizey = L{i}.blocks{1}.sizey + L{i}.blocks{2}.sizey;
                        new_block.size = sqrt(new_block.sizex) * sqrt(new_block.sizey);
                        XX = horzcat(L{i}.blocks{1}.X, L{i}.blocks{2}.X);
                        YY = horzcat(L{i}.blocks{1}.Y, L{i}.blocks{2}.Y);
                        if size(XX, 2) <= l
                            new_block.X = XX;
                            new_block.Y = YY;
                        else
                            [new_block.X, new_block.Y] = cod(2, XX, YY, l); 
                        end
                        L{i}.blocks(1:2) = [];
                        L{i + 1}.blocks{end+1} = new_block;
                        
                    end
                end
            end
        end

        updateTime = updateTime + toc; % update end time

        % Query
        if t >= N && mod(t, gap) == 0
            Xw = X(:, max(1, t-N+1):t);
            Yw = Y(:, max(1, t-N+1):t);
            q = q + 1;
            % [eW, sizeW] = eh_query(L, l, t, X, Y, N, mx, my); 
            
            tic; % query start time
            qb = Block(mx, my, l, 0);
            qb.X = zeros(mx, 0);
            qb.Y = zeros(my, 0);
            sizeW = 0;
            for i = 1:length(L)
                for j = 1:length(L{i}.blocks)
                    sizeW = sizeW + 2*size(L{i}.blocks{j}.X,2);
                    [A, B] = merge(qb, L{i}.blocks{j}, l);
                    qb.X = A;
                    qb.Y = B;
                end
            end
            queryTime = queryTime + toc; % query end time
            
            sizeW = sizeW + 2*size(active.X,2); % buffer size
            % err = correlation_error(XW,YW,A,B); 
            XFYF = norm(Xw, 'fro') * norm(Yw, 'fro');
            eW = norm_XYT_ABT(Xw', Yw', qb.X', qb.Y')/XFYF;
            maxsize = max(maxsize, sizeW);
            err_avg = err_avg + eW;
            err_max = max(err_max, eW);
            clear Xw Yw;
            
            byteStream = getByteStreamFromArray(L);
            L_mem = length(byteStream);
            byteStream = getByteStreamFromArray(active);
            active_mem = length(byteStream);

            currentMemory = (L_mem + active_mem) / 1024 / 1024; 
            if currentMemory > maxMemoryUsage
                maxMemoryUsage = currentMemory;
            end

        end

        waitbar(t / n, h, sprintf('complete %d/%d (%.2f%%)', t, n, (t / n) * 100));

    end
    err_avg = err_avg / q;
    avgUpdateTime = updateTime / N;
    avgQueryTime = queryTime / q;
    fprintf('EH-COD, l=%d, avg_error %f, max_error %f, maxsize %d, update %f, query %f, memory %f\n', l, err_avg, err_max, maxsize, avgUpdateTime, avgQueryTime, maxMemoryUsage);
    result = [l, maxsize, err_avg, err_max, maxMemoryUsage, avgUpdateTime, avgQueryTime];
    close(h)
end

function [A, B] = merge(block1, block2, l)
    X = horzcat(block1.X, block2.X);
    Y = horzcat(block1.Y, block2.Y);
    if size(X, 2) <= l
        A = X;
        B = Y;
    else
        [A, B] = cod(2, X, Y, l); 
    end
end
