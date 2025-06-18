function [l_avg, l_ans, err_max, cost] = sampling(X, Y, l, N, gap, flag_time)
    tic;

    [mx, n] = size(X);
    my = size(Y, 1);
    l_ans = -1;
    l_avg = l;
    err_avg = 0;
    err_max = 0;

    Q = {};
    digit = 0;

    for t = 1:n

        % Remove expiring blocks
        while ~isempty(Q) && Q{1}.t < t - N
            Q(1) =  [];
        end

        x = X(:, t);
        y = Y(:, t);
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
        
        % Query
        if t >= N && mod(t, gap) == 0 && flag_time == 0
            digit = digit + 1;
            err = query(t, Q, X, Y, N, l, mx, my);
            err_avg = err_avg + err;
            err_max = max(err_max,err);
            l_ans = max(l_ans, length(Q));
        end
    end

    err_avg = err_avg / digit;
    cost = toc;
end

function weight = calculate_weight(x, y)
    u = rand();
    x2 = norm(x, 2);
    y2 = norm(y, 2);
    weight = (x2 * y2) / u ;
end

function err = query(t, Q, X, Y, N, l, mx, my)
    XW = X(:, max(1, t-N+1):t);
    YW = Y(:, max(1, t-N+1):t);
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

    XFYF = XWF * YWF;
    err = norm_XYT_ABT(XW', YW', A', B') / XFYF;
end


