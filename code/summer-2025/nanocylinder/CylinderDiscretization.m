% Testing discretization

tic; % start timer to measure run time

[lambda12, sens12] = CylinderLambdaSens(1.33, 1.41, 12, [20, 20], 600, 3);
[lambda30, sens30] = CylinderLambdaSens(1.33, 1.41, 30, [20, 20], 600, 3);
[lambda40, sens40] = CylinderLambdaSens(1.33, 1.41, 40, [20, 20], 600, 3);
[lambda50, sens50] = CylinderLambdaSens(1.33, 1.41, 50, [20, 20], 600, 3);
[lambda60, sens60] = CylinderLambdaSens(1.33, 1.41, 60, [20, 20], 600, 3);

lambda = [lambda12, lambda30, lambda40, lambda50, lambda60];
sens = [sens12, sens30, sens40, sens50, sens60];
N = [lambda; sens];
M = transpose(N);
writematrix(M, 'Discretization Data.xlsx', 'Sheet', 1, 'Range', 'A25:B29');


toc; % stop timer