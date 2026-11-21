%  Martinsson extinction spectra of Ag triangles  
function [lambda, sens, fractArea, sensDrop] = triangleLambdaDrop(nSub, nBulk, length, height, x0, stepSize, coeff)
    % substrate data
    [lambdaSub, sensSub] = subTriangleSens(nSub, nBulk, length, height, x0, stepSize, coeff);
    % Martinsson effective data
    epsEff = effectiveDielectricTriangle(nSub, nBulk, length, height);
    nEff = sqrt(epsEff);
    [lambdaEff, sensEff] = triangleLambdaSens(nEff, length, height, x0, stepSize, coeff);
    % suspension data
    [~, sensSol] = triangleLambdaSens(nBulk, length, height, x0, stepSize, coeff);
    % initialize variables to return
    lambda = [lambdaSub, lambdaEff];
    sens = [sensSub, sensEff, sensSol];
    fractArea = (sqrt(3) * length / 4)/(sqrt(3) * length / 2 + 3 * height);

    subDrop = (sensSub - sensSol) * 100/sensSub;
    effDrop = (sensEff - sensSol) * 100/sensEff;
    sensDrop = [subDrop, effDrop];
end




