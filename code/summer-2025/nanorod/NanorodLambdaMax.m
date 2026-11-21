function lspr = NanorodLambdaMax(n, diameter, height, x0, phiSteps, thetaSteps, stepSize)
    mat = [Material(n ^ 2, 1), Material( epstable( 'gold.dat' ), 1 )];
    rules = quadboundary.rules( 'quad3', triquad( 1 ) );
    
    steps = ceil(height/stepSize); % discretization steps
    p = trirod(diameter, height, [phiSteps, thetaSteps, steps]);
    p.verts(:,3) = p.verts(:,3) - min(p.verts(:,3)) + 1e-2;
    tau = BoundaryEdge(mat, p, [1,2]);
    
    bem = galerkin.bemsolver(tau, 'rules', rules, 'waitbar', 1);

    % planewave excitation
    exc = galerkin.planewave( [ 0, 0, 1 ], [ 1, 0, 0 ] );

    [x, ~]=fminunc(@(x)-extinction(exc, bem \ exc(tau, 2 * pi ./ x)) ,x0); % find peak wavelength
    lspr = x;
end
