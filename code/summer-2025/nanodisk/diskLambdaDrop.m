%  Martinsson extinction spectra of Au triangles  
function [lambdaSub, lambda, sens, fractArea] = diskLambdaDrop(nSub, nBulk, diameter, height, x0, x0Sub, stepSize, vertices)
    % substrate data
    [lambdaSub, sensSub] = subDiskSens(nSub, nBulk, diameter, height, x0Sub, stepSize, vertices);
    % suspension data
    [lambda, sensSol] = diskLambdaSens(nBulk, diameter, height, x0, stepSize, vertices);
    % initialize variables to return
    sens = [sensSub, sensSol];
    radius = diameter/2;
    fractArea = radius/(2*radius + 2*height);
end
