function [err_avg, err_max, maxsize] = sampling_time(X, Y, timestamps, N, l, gap)
    start_time = tic;

    [mx, n] = size(X);
    my = size(Y, 1);
    maxsize = -1;
    err_avg = 0;
    err_max = 0;

    Q = {};
    digit = 0;
    
    t_max = timestamps(end);
    current_data_idx = 1;
    for t = 0:t_max

        % Remove expiring blocks
        while ~isempty(Q) && Q{1}.t < t - N
            Q(1) =  [];
        end
        
        while current_data_idx <= n && timestamps(current_data_idx) == t
            x = X(:, current_data_idx);
            y = Y(:, current_data_idx);
            current_data_idx = current_data_idx + 1;
            weight = calculate_weight(x, y);
            new_candidate = Candidate(x, y, t, weight);
            % Update queue and rank
            for i = 1:length(Q)
                if weight > Q{i}.weight
                    Q{i}.rank = Q{i}.rank + 1;
                end
            end
            % Retain elements ranked within 1
            tempQ = Q;
            Q = {};
            for i = 1:length(tempQ)
                if tempQ{i}.rank <= l
                    Q{end+1} = tempQ{i};
                end
            end
            Q{end+1} = new_candidate;
            
        end
        
        % Query
        if t >= N && mod(t, gap) == 0
            window_start_time = t - N + 1;
            window_indices = find(timestamps >= window_start_time & timestamps <= t);
            Xw = X(:, window_indices);
            Yw = Y(:, window_indices);
            digit = digit + 1;
            err = query(Q, Xw, Yw, l, mx, my);
            err_avg = err_avg + err;
            err_max = max(err_max,err);
            maxsize = max(maxsize, length(Q));
        end
    end

    err_avg = err_avg / digit;
    end_time = toc(start_time);
    fprintf('Sampling, for l=%d, the avg_error is %f, the max_error is %f, maxsize is %d, time cost %f\n', l, err_avg, err_max, maxsize, end_time);
end

function weight = calculate_weight(x, y)
    u = rand();
    x2 = norm(x, 2);
    y2 = norm(y, 2);
    weight = (x2 * y2) / u ;
end

function err = query(Q, XW, YW, l, mx, my)
    XWF = norm(XW, 'fro');
    YWF = norm(YW, 'fro');

    % Extract weights and sort them in descending order
    weights = arrayfun(@(c) c.weight, [Q{:}]);
    [~, idx] = sort(weights, 'descend');
    Q_sorted = Q(idx);

    % Select the element with the highest weight
    top_l_elements = Q_sorted(1:min(l, numel(Q_sorted)));

    % Initialise two matrices to store the x and y vectors
    A = zeros(mx, l);
    B = zeros(my, l);

    % Filling Matrix
    for i = 1:length(top_l_elements)
        candidate = top_l_elements{i}; 
        x = candidate.x;
        y = candidate.y;
        a = XWF / (sqrt(l) * norm(x, 2));
        b = YWF / (sqrt(l) * norm(y, 2));
        A(:, i) = a * x;  
        B(:, i) = b * y; 
    end

    err = correlation_error(XW,YW,A,B); 
end


