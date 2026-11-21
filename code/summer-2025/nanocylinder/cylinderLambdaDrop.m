%  Martinsson extinction spectra of Au triangles  
function [lambda, sens] = cylinderLambdaDrop(nSub, nBulk, diameter, height, x0, stepSize, vertices)
    % substrate data
    [lambdaSub, sensSub] = subCylinderSens(nSub, nBulk, diameter, height, x0, stepSize, vertices);
    % Martinsson effective data
    % epsEff = effectiveDielectricTriangle(nSub, nBulk, length, height);
    % nEff = sqrt(epsEff);
    % [lambdaEff, sensEff] = triangleLambdaSens(nEff, length, height, x0, stepSize, coeff);
    % suspension data
    [~, sensSol] = CylinderLambdaSens(nBulk, diameter, height, x0, stepSize, vertices);
    % initialize variables to return
    % lambda = [lambdaSub, lambdaEff];
    % sens = [sensSub, sensEff, sensSol];
    lambda = lambdaSub;
    sens = [sensSub, sensSol];
    

    % subDrop = (sensSub - sensSol) * 100/sensSub;
    % effDrop = (sensEff - sensSol) * 100/sensEff;
    % sensDrop = [subDrop, effDrop];
end




