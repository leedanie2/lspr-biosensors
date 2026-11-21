function plotTriangleSpectra(n, length, height, stepSize, coeff)
    mat = [Material(n ^ 2, 1), Material( epstable( 'silver.dat' ), 1 )];

    steps = ceil(height/stepSize); % discretization steps, each step is 2 nm
    edge = edgeprofile( height, steps, 'mode', '11');
    poly = polygon(3, 'size', [(sqrt(3)/2) * length, length]);
    refun = @(pos,d) coeff*ones(size(pos(:,1))); % refine mesh of triangle
    p = tripolygon(poly, edge, 'refun', refun); % extrude polygon to make nanotriangle
    p.verts( :, 3 ) = p.verts( :, 3 ) - min(p.verts( :, 3)) + 1e-2;

    exc = galerkin.planewave( [ 0, 1, 0 ], [ 0, 0, 1 ] );
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
    ylabel( 'Extinction Efficiency' );
end