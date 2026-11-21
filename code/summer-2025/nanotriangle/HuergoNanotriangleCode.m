% Huergo 2016 NT data

% Gold NTs in water
mat = [Material(1.33 ^ 2, 1), Material( epstable( 'gold.dat' ), 1 )];

% NT dimensions
sideLength = [40, 60, 80, 100, 120];
depth = 7;

rules = quadboundary.rules( 'quad3', triquad( 1 ) );
steps = ceil(depth/3) + 1;
tau = cell(1, 5);
bem = cell(1, 5);

for i = 1:5
    coeff = 2.5 * (sideLength(i)/40) ^ 2;
    refun = @(pos,d) coeff*ones(size(pos(:,1))); % refine mesh of triangle
    poly = polygon(3, 'size', [(sqrt(3)/2) * sideLength(i), sideLength(i)]);
    edge = edgeprofile(depth, steps, 'mode', '11');
    p = tripolygon(poly, edge, 'refun', refun);

    tau{i} = BoundaryEdge(mat, p, [1,2]);
    bem{i} = galerkin.bemsolver(tau{i}, 'rules', rules, 'waitbar', 1);
end
% planewave excitation
exc = galerkin.planewave( [ 0, 1, 0 ], [ 0, 0, 1 ] );

%  light wavelength in vacuum
lambda = linspace( 400, 1200, 100 );
k0 = 2 * pi ./ lambda;

nk0 = numel(k0);
npol = size(exc.pol, 1);
cext = zeros(nk0, npol, 5);

for j = 1 : nk0
    for k = 1:5
    %  solution of BEM equations
    sol = bem{k} \ exc( tau{k}, k0( j ) );

    %  extinction cross sections
    cext( j, : , k) = extinction( exc, sol );
    end
  multiWaitbar( 'BEM solver', j / numel( k0 ) );
end

%  close waitbar
multiWaitbar( 'CloseAll' );

%% final plot of figure 4.b)
figure(1); hold on
plot(lambda, cext(:, 1, 1), 'g--'); 
plot(lambda, cext(:, 1, 2), 'g-');   
plot(lambda, cext(:, 1, 3), 'y-');  
plot(lambda, cext(:, 1, 4), '-',  'Color', [1, 0.533, 0]);
plot(lambda, cext(:, 1, 5), 'r-');  

hold off