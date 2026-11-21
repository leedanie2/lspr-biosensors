function plotCylinderSpectra(n, diameter, height, stepSize, vertices)
    mat = [Material(n ^ 2, 1), Material( epstable( 'gold.dat' ), 1 )];

    steps = ceil(height/stepSize); % discretization steps
    edge = edgeprofile( height, steps, 'mode', '11');
    poly = polygon(vertices, 'size', [diameter, diameter]);
    p = tripolygon(poly, edge);

    exc = galerkin.planewave( [ 0, 0, 1 ], [ 1, 0, 0 ] );
    tau = BoundaryEdge(mat, p, [1,2]);
    rules = quadboundary.rules( 'quad3', triquad( 1 ) );
    bem = galerkin.bemsolver(tau, 'rules', rules, 'waitbar', 1);
    %  light wavelength in vacuum
    lambda = linspace( 400, 900, 50 );
    k0 = 2 * pi ./ lambda;

    % Initialize before looping over wavenumbers
    nk0 = numel(k0);
    npol = size(exc.pol, 1);
    cext = zeros(nk0, npol );
    %% 

    %  loop over wavenumbers
    for i = 1 : nk0
      %  solution of BEM equations
      sol = bem \ exc( tau, k0( i ) );
      %  extinction cross sections
      cext( i, : ) = extinction( exc, sol );
    end
    %% 
    
    %%  final plot for Figure 1
    figure(1);
    plot(lambda, cext,'b'); hold off
    xlabel( 'Wavelength, nm' );
    ylabel( 'Extinction' );
end