% disk in suspension

function lspr = diskLambdaMax(n, diameter, height, x0, stepSize, vertices)
    mat = [Material(n ^ 2, 1), Material( epstable( 'gold.dat' ), 1 )];
    rules = quadboundary.rules( 'quad3', triquad( 1 ) );

    steps = ceil(height/stepSize); % discretization steps
    edge = edgeprofile( height, steps, 'mode', '11');
    poly = polygon(vertices, 'size', [diameter, diameter]);
    p = tripolygon(poly, edge);
    p.verts( :, 3 ) = p.verts( :, 3 ) - min(p.verts( :, 3)) + 1e-2;

    tau = BoundaryEdge(mat, p, [1,2]);
    bem = galerkin.bemsolver(tau, 'rules', rules, 'waitbar', 1);

    %  planewave excitation
    exc = galerkin.planewave( [ 0, 1, 0 ], [ 0, 0, 1 ] );
    options = optimset('Display', 'iter', 'MaxIter', 100, 'TolX', 1e-6);
    [lspr, ~]=fminsearch(@(x)-extinction(exc, bem \ exc(tau, 2 * pi ./ x)) ,x0, options); % find peak wavelength
end