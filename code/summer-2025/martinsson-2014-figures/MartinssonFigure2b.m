% Martinsson Figure 2(b)/2(e) Ag Nanocube Sensitivity in solution vs substrate

 %  nanocube
length = 40;
poly = polygon(4, 'size', [length, length]);
stepSize = 5;
steps = ceil(length/stepSize) + 1;
edge = edgeprofile(length, steps, 'mode', '11');
refun = @(pos,d) 3*ones(size(pos(:,1)));
p = tripolygon(poly, edge, 'refun', refun);
p.verts( :, 3 ) = p.verts( :, 3 ) - min(p.verts( :, 3)) + 1e-2;
%  material properties of layer structure
mat1 = Material( 1.52 ^ 2, 1 ); % substrate
mat0 = Material( 1.33 ^ 2, 1); % relative baseline
mat3 = Material( epstable( 'silver.dat' ), 1 );
rules = quadboundary.rules( 'quad3', triquad( 1 ) );

n = zeros(1,10);
deltaLambdaSub = zeros(1,10);
deltaLambdaSol = zeros(1,10);

index = 1;
for i = 1.33:0.01:1.42
    n(index) = i;
    [ deltaLambdaSub(index), deltaLambdaSol(index) ] = cubeLambdaDrop(p, i, 450, mat1, mat0, mat3, rules);
    index = index + 1;
end

%% 

plot(n, deltaLambdaSub, 'r'); hold on
plot(n, deltaLambdaSol, 'k'); hold off

xlim([1.33 1.42]);
ylim([0 30]);
xlabel( 'Refractive Index' );
ylabel( 'Δλmax / nm' );