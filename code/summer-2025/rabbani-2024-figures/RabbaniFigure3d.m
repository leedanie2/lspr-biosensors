% nanotriangle (40 nm x 10 nm)
sideLength = 40;
depth = 10;
poly = polygon(3);
poly = scale(poly, sideLength/sqrt(3));
% analyte layer 10 nm
sideLength10 = sideLength + sideLength/sqrt(3);
depth10 = depth + 20;
poly10 = polygon(3);
poly10 = scale(poly10, sideLength10/sqrt(3));
% analyte layer 20 nm
sideLength20 = sideLength + 2 * sideLength/sqrt(3);
depth20 = depth + 40;
poly20 = polygon(3);
poly20 = scale(poly20, sideLength20/sqrt(3));

% toggle disc steps, # lambda points, vertices to change accuracy/runtime
steps = ceil(depth/5) + 1;
edge = edgeprofile(depth, steps, 'mode', '11');
refun = @(pos,d) 5*ones(size(pos(:,1))); % refine mesh of triangle
p1 = tripolygon(poly, edge, 'refun', refun); % extrude polygon to make nanotriangle

steps10 = ceil(depth10/5) + 1;
edge10 = edgeprofile(depth10, steps10, 'mode', '11');
p2 = tripolygon(poly10, edge10, 'refun', refun);

steps20 = ceil(depth20/5) + 1;
edge20 = edgeprofile(depth20, steps20, 'mode', '11');
p3 = tripolygon(poly20, edge20, 'refun', refun); 

% dielectric functions (outward to inward)
mat1 = Material( 1.33^2, 1);
mat2a = Material( 1.40^2, 1);
mat2b = Material( 1.41^2, 1);
mat3 = Material(epstable('gold.dat'), 1);

mat_a = [mat1, mat2a, mat3];
mat_b = [mat1, mat2b, mat3];

% boundary elements
tau = BoundaryEdge(mat_a, p1, [1,3]); % no analyte
tau_a140 = BoundaryEdge(mat_a, p2, [2,1], p1, [3, 2]); % 10 nm analyte, RI = 1.40
tau_a141 = BoundaryEdge(mat_b, p2, [2,1], p1, [3,2]); % 10 nm analyte, RI = 1.41
tau_b140 = BoundaryEdge(mat_a, p3, [2,1], p1, [3,2]); % 20 nm analyte, RI = 1.40
tau_b141 = BoundaryEdge(mat_b, p3, [2,1], p1, [3, 2]); % 20 nm analyte, RI = 1.41

% initialize bem solver
rules = quadboundary.rules( 'quad3', triquad( 1 ) );
bem = galerkin.bemsolver(tau, 'rules', rules, 'waitbar', 1);
bem_a140 = galerkin.bemsolver(tau_a140, 'rules', rules, 'waitbar', 1);
bem_a141 = galerkin.bemsolver(tau_a141, 'rules', rules, 'waitbar', 1);
bem_b140 = galerkin.bemsolver(tau_b140, 'rules', rules, 'waitbar', 1);
bem_b141 = galerkin.bemsolver(tau_b141, 'rules', rules, 'waitbar', 1);

% planewave excitation
exc = galerkin.planewave( [ 0, 1, 0 ], [ 0, 0, 1 ] );

%  light wavelength in vacuum
lambda = linspace( 650, 750, 30);
k0 = 2 * pi ./ lambda;

%  allocate optical cross sections
[cext, cext_a140, cext_a141, cext_b140, cext_b141] = deal( zeros( numel( k0 ), size( exc.pol, 1 ) ) );
multiWaitbar( 'BEM solver', 0, 'Color', 'g', 'CanCancel', 'on' );

%  loop over wavenumbers
for i = 1 : numel( k0 )
  %  solution of BEM equations
  sol = bem \ exc( tau, k0( i ) );
  sol_a140 = bem_a140 \ exc( tau_a140, k0( i ) );
  sol_a141 = bem_a141 \ exc( tau_a141, k0( i ) );
  sol_b140 = bem_b140 \ exc( tau_b140, k0( i ) );
  sol_b141 = bem_b141 \ exc( tau_b141, k0( i ) );

  %  extinction cross sections
  cext( i, : ) = extinction( exc, sol );
  cext_a140( i, : ) = extinction( exc, sol_a140 );
  cext_a141( i, : ) = extinction( exc, sol_a141 );
  cext_b140( i, : ) = extinction( exc, sol_b140 );
  cext_b141( i, : ) = extinction( exc, sol_b141 );
 
  multiWaitbar( 'BEM solver', i / numel( k0 ) );
end
%  close waitbar
multiWaitbar( 'CloseAll' );

% geometric cross section
cgeo = pi * (3000 * sqrt(3)/pi)^(2/3);
% extinction efficiencies
exteff = cext/cgeo;
exteff_a140 = cext_a140/cgeo;
exteff_a141 = cext_a141/cgeo;
exteff_b140 = cext_b140/cgeo;
exteff_b141 = cext_b141/cgeo;

%%  final plot
figure(1); 
plot( lambda, exteff, 'k'); hold on
plot( lambda, exteff_a140, 'r');  hold on
plot( lambda, exteff_a141, 'b');  hold on
plot( lambda, exteff_b140, 'r--');  hold on
plot( lambda, exteff_b141, 'b--');  hold off

xlabel( 'Wavelength (nm)' );
ylabel( 'Extinction Efficiency' );
