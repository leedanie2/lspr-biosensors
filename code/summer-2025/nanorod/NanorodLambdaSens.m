function [lambda, sens] = NanorodLambdaSens(n1, n2, diameter, height, x0, phiSteps, thetaSteps, stepSize)
    lambdaMax2 = NanorodLambdaMax(n2, diameter, height, x0, phiSteps, thetaSteps, stepSize);
    lambda = NanorodLambdaMax(n1, diameter, height, x0, phiSteps, thetaSteps, stepSize);
    sens = (lambdaMax2 - lambda)/(n2 - n1);
end