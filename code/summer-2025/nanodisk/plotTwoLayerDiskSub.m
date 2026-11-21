function plotTwoLayerDiskSub(nBulk, nSub, nLocal, diameter, height, analyteThickness, stepSize, vertices, eps)
   
    % Define outer and inner radii
    diameterInner = diameter + eps;
    diameterOuter = diameter + 2 * analyteThickness;
    
    % Edge profile for particle
    steps = ceil(height / stepSize) + 1;
    edge = edgeprofile(height, steps, 'mode', '11');

    % Edge profile for the cap (top disk)
    heightCap = analyteThickness;
    stepSizeAn = stepSize + 3;
    stepsCap = ceil(heightCap / stepSizeAn) + 1;
    edgeCap = edgeprofile(heightCap, stepsCap, 'mode', '11');

    % Make ring polygon (bottom half of analyte layer)
    polyRingOuter = polygon(vertices, 'size', [diameterOuter, diameterOuter], 'dir', +1);
    polyRingInner = polygon(vertices, 'size', [diameterInner, diameterInner], 'dir', -1);
    polyRing = [polyRingOuter, polyRingInner];
    edgeRing = edgeprofile(height, stepSizeAn, 'mode', '11');
    pRing = tripolygon(polyRing, edgeRing); 

    % Create the cap (top half of analyte layer)
    polyCap = polygon(vertices, 'size', [diameterOuter, diameterOuter], 'dir', +1);
    pCap = tripolygon(polyCap, edgeCap);

    % Shift the cap to be on top of the ring
    ringTop = max(pRing.verts(:,3));
    capBottom = min(pCap.verts(:,3));
    pCap = shift(pCap, [0, 0, ringTop - capBottom + eps]);
    
    % Create the inner metallic disk (p)
    poly = polygon(vertices, 'size', [diameter, diameter]);
    p = tripolygon(poly, edge); % Using 'edge' for height of 20, same as original disk height

    % Shift the metallic disk relative to the analyte layer (bottom of pRing)
    ringBottom = min(pRing.verts(:,3));
    npBottom = min(p.verts(:,3));
    shiftRing = (ringBottom - npBottom);
    p  = shift(p, [0, 0, shiftRing]);

    % Shift everything so the bottom is right above substrate
    lowest = min([min(pRing.verts(:,3)), min(pCap.verts(:,3)), min(p.verts(:,3))]);
    vertShift = max(0, -lowest + 1e-2);
    pRing = shift(pRing, [0, 0, vertShift]);
    pCap = shift(pCap, [0, 0, vertShift]);
    p = shift(p, [0, 0, vertShift]);

    % Define materials
    mat1 = Material( nSub ^ 2, 1 ); % Substrate
    mat2 = Material( nBulk ^ 2, 1 ); % Bulk
    mat3 = Material( epstable( 'gold.dat' ), 1 ); 
    mat4a = Material( nLocal ^ 2, 1 ); % Analyte
    
    matA = [mat1, mat2, mat3, mat4a];
    
    % Define stratified layer structure (substrate and bulk)
    layer = stratified.layerstructure( matA( 1 : 2 ), 0 ); % Substrate (mat1) and Bulk (mat2)

    tau = BoundaryEdge(matA, pRing, [4, 2], pCap, [4, 2], p, [3, 4]); 
    % plot(tau, 'EdgeColor', 'k');
    bem = stratified.bemsolver( tau, layer, 'waitbar', 1 );

    %  planewave excitation
    exc = stratified.planewave( layer, [ 0, 1, 0 ], [ 0, 0, -1 ] );

    %  light wavelength in vacuum
    lambda = linspace( 400, 800, 20 );
    k0 = 2 * pi ./ lambda;

    % Initialize before looping over wavenumbers
    nk0 = numel(k0);
    npol = size(exc.pol, 1);
    cext = zeros(nk0, npol );

    %  loop over wavenumbers
    for i = 1 : nk0
      %  solution of BEM equations
      sol = bem \ exc( tau, k0( i ) );
      %  extinction cross sections
      cext( i, : ) = extinction( exc, sol );
    end

    %%  final plot
    figure();
    plot(lambda, cext,'b'); hold off
    xlabel( 'Wavelength, nm' );
    ylabel( 'Extinction' );
end
