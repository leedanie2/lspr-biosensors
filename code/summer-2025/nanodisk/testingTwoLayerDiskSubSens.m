function sens = testingTwoLayerDiskSubSens(nBulk, nSub, nLocal, diameter, height, analyteThickness, x0, stepSize, vertices, eps)

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
    anVerts = vertices - 4;
    polyRingOuter = polygon(anVerts, 'size', [diameterOuter, diameterOuter], 'dir', +1);
    polyRingInner = polygon(vertices, 'size', [diameterInner, diameterInner], 'dir', -1);
    polyRing = [polyRingOuter, polyRingInner];
    edgeRing = edgeprofile(height, stepSizeAn, 'mode', '11');
    pRing = tripolygon(polyRing, edgeRing); 

    % Create the cap (top half of analyte layer)
    polyCap = polygon(anVerts, 'size', [diameterOuter, diameterOuter], 'dir', +1);
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

    % Shift everything so the bottom is 0.01 nm above substrate
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

    % Define incident plane wave
    exc = stratified.planewave( layer, [ 0, 1, 0 ], [ 0, 0, -1 ] );
    
    % Optimization options for fminsearch
    options = optimset('Display', 'iter', 'MaxIter', 50, 'TolX', 1e-3);
    tau = BoundaryEdge(matA, pRing, [4, 2], pCap, [4, 2], p, [3, 4]); 
    bem = stratified.bemsolver( tau, layer, 'waitbar', 1 );
    [lambda1, ~] = fminsearch(@(x) -extinction(exc, bem \ exc(tau, 2 * pi ./ x)), x0, options);
    
    % finding lspr wavelength 2
    nBulk2 = nBulk + 0.01;
    mat2b = Material( nBulk2 ^ 2, 1);
    nLocal2 = nLocal + 0.01;
    mat4b = Material( nLocal2 ^ 2, 1 );
    matB = [mat1, mat2b, mat3, mat4b];
   
    tau = BoundaryEdge(matB, pRing, [4, 2], pCap, [4, 2], p, [3, 4]);
    bem = stratified.bemsolver( tau, layer, 'waitbar', 1 );
    [lambda2, ~] = fminsearch(@(x) -extinction(exc, bem \ exc(tau, 2 * pi ./ x)), lambda1, options);
    
    % Calculate sensitivity
    sens = (lambda2 - lambda1) / 0.01;

end