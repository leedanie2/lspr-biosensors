function [lambda, sens] = triangleLambdaSens(n, length, height, x0, stepSize, coeff)
    n2 = n + 0.01;
    lambdaMax2 = triangleLambdaMax(n2, length, height, x0, stepSize, coeff);
    lambda = triangleLambdaMax(n, length, height, x0, stepSize, coeff);
    sens = (lambdaMax2 - lambda)/(0.01);
end