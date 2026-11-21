function [lambda, sens] = searchLambdaSens(n1, n2, height, dimensions, x0)
    lambdaMax2 = searchLambdaMax(n2, height, dimensions, x0);
    lambda = searchLambdaMax(n1, height, dimensions, x0);
    sens = (lambdaMax2 - lambda)/(n2 - n1);
end