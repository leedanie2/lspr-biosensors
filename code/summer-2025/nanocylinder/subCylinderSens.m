function [lambda, sens] = subCylinderSens(nSub, nBulk, diameter, height, x0, stepSize, vertices)
    nBulk2 = nBulk + 0.01;
    lambdaMax2 = subCylinderLambdaMax(nSub, nBulk2, diameter, height, x0, stepSize, vertices);
    lambda = subCylinderLambdaMax(nSub, nBulk, diameter, height, x0, stepSize, vertices);
    
    sens = (lambdaMax2 - lambda)/(0.01);
end