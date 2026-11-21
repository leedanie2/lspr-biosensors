function lspr = searchLambdaMax(n, height, dimensions, x0)
    mat = [Material(n ^ 2, 1), Material( epstable( 'gold.dat' ), 1 )];
    rules = quadboundary.rules( 'quad3', triquad( 1 ) );
    % steps = ceil(height/stepSize);
    edge = edgeprofile( height, 9, 'mode', '11');
    poly = polygon(10, 'size', dimensions);
    p = tripolygon( poly, edge);
    p.verts(:,3) = p.verts(:,3) - min(p.verts(:,3)) + 1e-2;
    tau = BoundaryEdge(mat, p, [1,2]);
    bem = galerkin.bemsolver(tau, 'rules', rules, 'waitbar', 1);

    % planewave excitation
    exc = galerkin.planewave( [ 0, 0, 1 ], [ 1, 0, 0 ] );

    options = optimset('Display', 'iter', 'MaxIter', 100, 'TolX', 1e-4);
    [x, ~]=fminsearch(@(x)-extinction(exc, bem \ exc(tau, 2 * pi ./ x)) ,x0, options); % find peak wavelength
    lspr = x;
end