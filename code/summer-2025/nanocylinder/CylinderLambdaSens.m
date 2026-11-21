% Sensitivity and lambda max of n1 for gold nanocylinders

function [lambda, sens] = CylinderLambdaSens(n, diameter, height, x0, stepSize, vertices)
    n2 = n + 0.01;
    lambdaMax2 = CylinderLambdaMax(n2, diameter, height, x0, stepSize, vertices);
    lambda = CylinderLambdaMax(n, diameter, height, x0, stepSize, vertices);
    sens = (lambdaMax2 - lambda)/0.01;
end