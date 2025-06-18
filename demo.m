clear;
warning('off', 'all');

dataset = 'syn_uni.mat';
load(dataset);

[L_estimate, R_estimate] = checkR(X, Y);
fprintf('L_estimate = %f R_estimate = %f\n', L_estimate, R_estimate);

% Set window size N and query gap
N = 4000; 
gap = 400;

% Set sketch size
ls = [10,25,40];

[mx, my, n, R] = info(X, Y);
fprintf('X size %d × %d, Y size %d × %d, N = %d, R = %f \n', mx, n, my, n, N, R);

logTimeStr = datestr(now, 'yyyymmdd_HHMM');
[~, dataset_name, ~] = fileparts(dataset);  

logFile = sprintf('output/log_%s_%s.txt', dataset_name, logTimeStr);
fid = fopen(logFile, 'a'); 

fprintf(fid, 'Dataset: %s | X: %d×%d, Y: %d×%d, N = %d, R = %.6f, L_estimate = %.6f, R_estimate = %.6f\n', ...
    dataset, mx, n, my, n, N, R, L_estimate, R_estimate);

% hDS-COD
for k = 1:length(ls)
    result = hds_cod(X, Y, N, ls(k), L_estimate, R_estimate, gap);
    fprintf(fid, 'hDS-COD, l=%d, maxsize %.6f, avg_error %.6f, max_error %.6f, memory %.6f, update %.6f, query %.6f\n', ...
        result(1), result(2), result(3), result(4), result(5), result(6), result(7));
end

% aDS-COD
for k = 1:length(ls)
    result = ads_cod(X, Y, N, ls(k), R_estimate, gap);
    fprintf(fid, 'aDS-COD, l=%d, maxsize %.6f, avg_error %.6f, max_error %.6f, memory %.6f, update %.6f, query %.6f\n', ...
        result(1), result(2), result(3), result(4), result(5), result(6), result(7));
end

% EH-COD
for k = 1:length(ls)
    result = eh_cod(X, Y, N, ls(k), gap);
    fprintf(fid, 'EH-COD, l=%d, maxsize %.6f, avg_error %.6f, max_error %.6f, memory %.6f, update %.6f, query %.6f\n', ...
        result(1), result(2), result(3), result(4), result(5), result(6), result(7));
end

% DI-COD
L = floor(log2(R_estimate));
for i = L-4:L+3
    try
        result = di_cod(X, Y, N, i, gap);
        fprintf(fid, 'DI-COD, level=%d, maxsize %.6f, avg_error %.6f, max_error %.6f, memory %.6f, update %.6f, query %.6f\n', ...
            result(1), result(2), result(3), result(4), result(5), result(6), result(7));
    catch ME
    end
end

fclose(fid);


