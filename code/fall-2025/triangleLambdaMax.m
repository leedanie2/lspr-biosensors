function lspr = triangleLambdaMax(n, length, height, x0, stepSize, coeff)
    mat = [Material(n ^ 2, 1), Material( epstable( 'gold.dat' ), 1 )];
    rules = quadboundary.rules( 'quad3', triquad( 1 ) );
    steps = ceil(height/stepSize); 

    % equilateral triangle
    poly = polygon(3, 'size', [(sqrt(3)/2) * length, length]);
    edge = edgeprofile(height, steps, 'mode', '11');
    refun = @(pos,d) coeff*ones(size(pos(:,1))); % refine mesh of triangle
    p = tripolygon(poly, edge, 'refun', refun); % extrude polygon to make nanotriangle
    p.verts(:,3) = p.verts(:,3) - min(p.verts(:,3)) + 1e-2;

    tau = BoundaryEdge(mat, p, [1,2]);
    bem = galerkin.bemsolver(tau, 'rules', rules, 'waitbar', 1);
    % plot(tau, 'EdgeColor', 'b');

    % planewave excitation
    exc = galerkin.planewave( [ 0, 1, 0 ], [ 0, 0, 1 ] );


    options = optimset('Display', 'iter', 'MaxIter', 100, 'TolX', 1e-6);
    [lspr, ~]=fminsearch(@(x)-extinction(exc, bem \ exc(tau, 2 * pi ./ x)) ,x0, options); % find peak wavelength
end