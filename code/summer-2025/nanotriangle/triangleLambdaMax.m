function lspr = triangleLambdaMax(n, length, height, x0, stepSize, coeff)
    mat = [Material(n ^ 2, 1), Material( epstable( 'silver.dat' ), 1 )];
    rules = quadboundary.rules( 'quad3', triquad( 1 ) );
    steps = ceil(height/stepSize) + 1; 

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


    [x, ~]=fminunc(@(x)-extinction(exc, bem \ exc(tau, 2 * pi ./ x)) ,x0); % find peak wavelength
    lspr = x;
end