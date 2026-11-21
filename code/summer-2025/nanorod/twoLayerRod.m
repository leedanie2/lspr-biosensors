% sensitivity for Au nanorods in two layer systems
function [lambda_rod, sens] = twoLayerRod(nBulk, nLocal, t, height, diameter, x0, stepSize, anSteps)
    nLocal2 = nLocal + 0.01;
    % dielectric functions (outward to inward)
    mat1 = [Material(nBulk ^ 2, 1), Material(nLocal ^ 2, 1), Material( epstable( 'gold.dat' ), 1 )];
    mat2 = [Material(nBulk ^ 2, 1), Material(nLocal2 ^ 2, 1), Material( epstable( 'gold.dat' ), 1 )];
    % nanorod
    steps1 = ceil(height/stepSize) + 1;
    p1 = trirod(diameter, height, [16, 8, steps1]); 

    % analyte layer
    height_analyte = height + 2 * t;
    diameter_analyte = diameter + 2 * t;
    steps2 = ceil(height_analyte/anSteps) + 1;
    p2 = trirod(diameter_analyte, height_analyte, [16, 8, steps2]);

    % boundary element
    tau = BoundaryEdge(mat1, p2, [2,1], p1, [3, 2]);
    % initialize bem solver
    rules = quadboundary.rules( 'quad3', triquad( 1 ) );
    bem = galerkin.bemsolver(tau, 'rules', rules, 'waitbar', 1);
    % planewave excitation
    exc = galerkin.planewave( [ 0, 0, 1 ], [ 1, 0, 0 ] );
    % find peak wavelength
    [x1, ~]=fminunc(@(x)-extinction(exc, bem \ exc(tau, 2 * pi ./ x)) ,x0);
    lambda_rod = x1;
    % second LSPR calculation for sensitivity
    tau = BoundaryEdge(mat2, p2, [2,1], p1, [3, 2]);
    bem = galerkin.bemsolver(tau, 'rules', rules, 'waitbar', 1);
    [x2, ~]=fminunc(@(x)-extinction(exc, bem \ exc(tau, 2 * pi ./ x)) ,x0);
    % sensitivity
    sens = (x2 - x1)/(nLocal2 - nLocal);