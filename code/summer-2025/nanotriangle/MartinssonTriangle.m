% Martinsson effective medium approximation for nanotriangles

calcDrop = zeros(1,4);
epsEff = zeros(1,4);
effDrop = zeros(1,4);
%  material properties of layer structure
mat1 = Material( 1.5 ^ 2, 1 ); % glass substrate
mat2 = Material( 1.33 ^ 2, 1 ); % bulk medium (water)
mat3 = Material( epstable( 'silver.dat' ), 1 ); % silver NPs
%  set up layer structure
mat_sub = [ mat1, mat2, mat3 ];
mat_sol = [ mat2, mat3];
layer = stratified.layerstructure( mat_sub( 1 : 2 ), 0 );

% nanotriangles
lengths = [30, 40, 50, 60];
rules = quadboundary.rules( 'quad3', triquad( 1 ) );
stepSize = 3;
coeff = 2.5;
for i = 1:4
    epsEff(i) = effectiveDielectricTriangle(lengths(i), 10, 1.5, 1.33);
    [calcDrop(i), effDrop(i)] = triangleLambdaDrop(1.5, 1.33, epsEff(i), lengths(i), 10, 650, 3);
end
%% 

figure(1);
plot(lengths, calcDrop, 'k'); hold on
plot(lengths, effDrop, 'r'); hold off

xlabel( 'Triangle side length (nm)' );
ylabel( 'RI sensitivity drop proportion' );


