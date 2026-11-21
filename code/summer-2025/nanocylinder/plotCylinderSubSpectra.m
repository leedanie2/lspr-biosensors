function plotCylinderSubSpectra(nSub, nBulk, diameter, height, stepSize, vertices)
    mat1 = Material( nSub ^ 2, 1 ); % substrate
    mat2 = Material( nBulk ^ 2, 1 ); % bulk medium
    mat3 = Material( epstable( 'gold.dat' ), 1 ); % nanoparticle

    mat = [mat1, mat2, mat3];
    layer = stratified.layerstructure( mat( 1 : 2 ), 0 );

    steps = ceil(height/stepSize); % discretization steps
    edge = edgeprofile( height, steps, 'mode', '11');
    poly = polygon(vertices, 'size', [diameter, diameter]);
    p = tripolygon(poly, edge);

    % Rotate cylinder to lie flat along x-axis (long axis in x)
    theta = pi/2;  % 90 degree rotation about y-axis
    Ry = [cos(theta), 0, sin(theta); 0, 1, 0; -sin(theta), 0, cos(theta)];
    p.verts = (Ry * p.verts')';

    p.verts( :, 3 ) = p.verts( :, 3 ) - min(p.verts( :, 3)) + 1e-2;

    exc = stratified.planewave( layer, [ 1, 0, 0 ], [ 0, 0, 1 ] );
    tau = BoundaryEdge(mat, p, [3,2]);

    bem = stratified.bemsolver( tau, layer, 'waitbar', 1 ); %  initialize BEM solver
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