% sensitivity for Au nanocylinders in two layer systems
function [lambda_cyl, sens] = twoLayerCylinder(nBulk, nLocal, t, height, diameter, x0, stepSize, anSteps)
    nLocal2 = nLocal + 0.01;
    % dielectric functions (outward to inward)
    mat1 = [Material(nBulk ^ 2, 1), Material(nLocal ^ 2, 1), Material( epstable( 'gold.dat' ), 1 )];
    mat2 = [Material(nBulk ^ 2, 1), Material(nLocal2 ^ 2, 1), Material( epstable( 'gold.dat' ), 1 )];
    % nanocylinder
    steps1 = ceil(height/stepSize) + 1;
    edge1 = edgeprofile(height, steps1, 'mode', '11');
    poly1 = polygon(18, 'size', [diameter, diameter]); % 18 vertices
    p1 = tripolygon(poly1, edge1);
    % analyte layer
    height_analyte = height + 2 * t;
    diameter_analyte = diameter + 2 * t;
    % stepSize2 = stepSize + 2;
    steps2 = ceil(height_analyte/anSteps) + 1;
    edge2 = edgeprofile(height_analyte, steps2, 'mode', '11');
    poly2 = polygon(10, 'size', [diameter_analyte, diameter_analyte]);
    p2 = tripolygon(poly2, edge2);
    % boundary element
    tau = BoundaryEdge(mat1, p2, [2,1], p1, [3, 2]);
    % initialize bem solver
    rules = quadboundary.rules( 'quad3', triquad( 1 ) );
    bem = galerkin.bemsolver(tau, 'rules', rules, 'waitbar', 1);
    % planewave excitation
    exc = galerkin.planewave( [ 0, 0, 1 ], [ 1, 0, 0 ] );
    % find peak wavelength
    [x1, ~]=fminunc(@(x)-extinction(exc, bem \ exc(tau, 2 * pi ./ x)) ,x0);
    lambda_cyl = x1;
    % second peak wavelength for sensitivity calculations
    tau = BoundaryEdge(mat2, p2, [2,1], p1, [3, 2]);
    bem = galerkin.bemsolver(tau, 'rules', rules, 'waitbar', 1);
    [x2, ~]=fminunc(@(x)-extinction(exc, bem \ exc(tau, 2 * pi ./ x)) ,x0);
    % sensitivity
    sens = (x2 - x1)/(nLocal2 - nLocal);
