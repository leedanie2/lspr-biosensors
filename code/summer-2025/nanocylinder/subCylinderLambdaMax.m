% cylinder on substrate

function lspr = subCylinderLambdaMax(nSub, nBulk, diameter, height, x0, stepSize, vertices)
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

    options = optimset('Display', 'iter', 'MaxIter', 100, 'TolX', 1e-2);
    [lspr, ~]=fminsearch(@(x)-extinction(exc, bem \ exc(tau, 2 * pi ./ x)) ,x0, options); % find peak wavelength
end