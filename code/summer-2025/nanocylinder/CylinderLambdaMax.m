% Optimizing Lambda Max for gold nanocylinders
% n is refractive index of surrounding medium
% height is for height of cylinder
% dimensions is a 2x1 array for 2d cross-section (radius)
% x0 is initial guess for lambda max



function lspr = CylinderLambdaMax(n, diameter, height, x0, stepSize, vertices)
    mat = [Material(n ^ 2, 1), Material( epstable( 'gold.dat' ), 1 )];
    rules = quadboundary.rules( 'quad3', triquad( 1 ) );

    steps = ceil(height/stepSize); % discretization steps, each step is 2 nm
    edge = edgeprofile( height, steps, 'mode', '11');
    poly = polygon(vertices, 'size', [diameter, diameter]);
    p = tripolygon(poly, edge);
    tau = BoundaryEdge(mat, p, [1,2]);
    bem = galerkin.bemsolver(tau, 'rules', rules, 'waitbar', 1);

    % planewave excitation
    exc = galerkin.planewave( [ 0, 0, 1 ], [ 1, 0, 0 ] );

    options = optimset('Display', 'iter', 'MaxIter', 100, 'TolX', 1e-2);
    [x, ~]=fminsearch(@(x)-extinction(exc, bem \ exc(tau, 2 * pi ./ x)) ,x0, options); % find peak wavelength
    lspr = x;
end
