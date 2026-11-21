% Testing discretization for gold nanorods

tic; % start timer to measure run time

[lambda15, sens15] = NanorodLambdaSens(1.40, 1.41, 15, 40, 650, 12, 6, 3);
[lambda20, sens20] = NanorodLambdaSens(1.40, 1.41, 20, 40, 650, 12, 6, 3);
[lambda25, sens25] = NanorodLambdaSens(1.40, 1.41, 25, 40, 650, 12, 6, 3);


% % For making excel sheet overnight
% lambda = [lambda12, lambda30, lambda40, lambda50, lambda60];
% sens = [sens12, sens30, sens40, sens50, sens60];
% N = [lambda; sens];
% M = transpose(N);
% writematrix(M, 'Discretization Data.xlsx', 'Sheet', 1, 'Range', 'B2:C7');


toc; % stop timer