% Sensitivity and lambda max of n1 for gold nanocylinders

function [lambda, sens] = diskLambdaSens(n, diameter, height, x0, stepSize, vertices)
    n2 = n + 0.01;
    lambdaMax2 = diskLambdaMax(n2, diameter, height, x0, stepSize, vertices);
    lambda = diskLambdaMax(n, diameter, height, x0, stepSize, vertices);
    sens = (lambdaMax2 - lambda)/0.01;
end