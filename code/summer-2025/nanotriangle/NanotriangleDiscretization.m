% Testing discretization for gold nanotriangles

tic;

[lambda10, sens10] = triangleLambdaSens(1.40, 1.41, 10, 40, 700, 5);
[lambda15, sens15] = triangleLambdaSens(1.40, 1.41, 15, 40, 700, 5);
[lambda20, sens20] = triangleLambdaSens(1.40, 1.41, 20, 40, 700, 5);
toc;