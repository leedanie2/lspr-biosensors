% sensitivity for Au nanotriangles in two layer systems
function [lambda_tri, sens] = twoLayerTriangle(nBulk, nLocal, t, sideLength, depth, x0, stepSize,anSteps)
    nLocal2 = nLocal + 0.01;
    % dielectric functions (outward to inward)
    mat1 = [Material(nBulk ^ 2, 1), Material(nLocal ^ 2, 1), Material( epstable( 'gold.dat' ), 1 )];
    mat2 = [Material(nBulk ^ 2, 1), Material(nLocal2 ^ 2, 1), Material( epstable( 'gold.dat' ), 1 )];
    % nanotriangle
    steps = ceil(depth/stepSize) + 1;
    coeff1 = 3;
    refun = @(pos,d) coeff1*ones(size(pos(:,1))); % refine mesh of triangle
    poly = polygon(3, 'size', [(sqrt(3)/2) * sideLength, sideLength]);
    edge = edgeprofile(depth, steps, 'mode', '11');
    p1 = tripolygon(poly, edge, 'refun', refun);
    % analyte layer
    sideLength_analyte = sideLength + 2 * t * sqrt(3);
    depth_analyte = depth + 2 * t;
    % stepSize2 = stepSize + 2;
    steps2 = ceil(depth_analyte/anSteps) + 1;
    poly2 = polygon(3, 'size', [(sqrt(3)/2) * sideLength_analyte, sideLength_analyte]);
    edge2 = edgeprofile(depth_analyte, steps2, 'mode', '11');
    coeff2 = coeff1 * (sideLength_analyte/sideLength) ^ 2;
    refun2 = @(pos,d) coeff2*ones(size(pos(:,1)));
    p2 = tripolygon(poly2, edge2, 'refun', refun2);
    % boundary element
    tau = BoundaryEdge(mat1, p2, [2,1], p1, [3, 2]);
    % plot(tau, 'EdgeColor', 'b');

    % initialize bem solver
    rules = quadboundary.rules( 'quad3', triquad( 1 ) );
    bem = galerkin.bemsolver(tau, 'rules', rules, 'waitbar', 1);
    % planewave excitation
    exc = galerkin.planewave( [ 0, 1, 0 ], [ 0, 0, 1 ] );
    % find peak wavelength
    [x1, ~]=fminunc(@(x)-extinction(exc, bem \ exc(tau, 2 * pi ./ x)) ,x0);
    lambda_tri = x1;
    % second LSPR calculation
    tau = BoundaryEdge(mat2, p2, [2,1], p1, [3, 2]);
    bem = galerkin.bemsolver(tau, 'rules', rules, 'waitbar', 1);
    [x2, ~]=fminunc(@(x)-extinction(exc, bem \ exc(tau, 2 * pi ./ x)) ,x0);
    % sensitivity
    sens = (x2 - x1)/(nLocal2 - nLocal);

