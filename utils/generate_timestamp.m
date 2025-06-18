load("apr.mat")
[mx, my, n, R] = info(X, Y);

% Set the parameter for the Poisson process λ
lambda = 2;

% Generate n-1 time intervals, ensuring each interval is at least 1
inter_arrival_times = poissrnd(lambda, 1, n-1) + 1;

% Initialize an array for timestamps
timestamps = zeros(1, n);
timestamps(1) = 0; 
timestamps(2:end) = cumsum(inter_arrival_times);
disp('Generated timestamps:');
disp(timestamps);

save('apr_t.mat', 'X', 'Y', 'timestamps');
