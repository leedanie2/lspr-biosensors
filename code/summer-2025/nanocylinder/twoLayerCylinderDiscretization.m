% testing discretization for two layer Au nanocylinders (40 nm x 20 nm)
% analyte layer n = 1.40, bulk medium of water
% analyte thicknesses of 10 nm and 20 nm

tic;

[lambda10, sens10] = twoLayerCylinder(1.33, 1.40, 10, 40, 20, 650, 5, 10);
[lambda20, sens20] = twoLayerCylinder(1.33, 1.40, 20, 40, 20, 650, 5, 10);

toc;