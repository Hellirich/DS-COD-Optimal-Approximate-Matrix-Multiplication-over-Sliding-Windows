mx = 1000; 
my = 2000;
n = 10000; 

%% uniform
X = rand(mx, n);
Y = rand(my, n);
save('datasets/syn_uni.mat', 'X', 'Y');

%% EXP
% lambda_X = 1; 
% X = exprnd(1/lambda_X, [mx, n]);
% Y = exprnd(1/lambda_X, [my, n]); 
% save('datasets/syn_exp.mat', 'X', 'Y');

%% GAM
% k = 2; 
% theta = 1;
% X = gamrnd(k, theta, [mx, n]);
% Y = gamrnd(k, theta, [my, n]);
% save('datasets/syn_gamma.mat', 'X', 'Y');

%% syn_non_stationary
% X = zeros(mx, n); 
% Y = zeros(my, n); 
% for i = 1:n
%     mean_X = i / n * 5;  
%     std_X = 1 + 0.05 * (i/n); 
%     mean_Y = i / n * 5;  
%     std_Y = 1 + 0.05 * (i/n);
%     X(:, i) = mean_X + std_X * randn(mx, 1); 
%     Y(:, i) = mean_Y + std_Y * randn(my, 1); 
% end
% X=0.1*X;
% Y=0.1*Y;
% save('datasets/syn_non_stationary.mat', 'X', 'Y');

