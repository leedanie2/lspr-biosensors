% nanorod (40 nm x 20 nm)
height = 40;
diameter = 20;
% analyte layer 10 nm
height10 = height + 20;
diameter10 = diameter + 20;
% analyte layer 20 nm
height20 = height + 40;
diameter20 = diameter + 40;

% toggle disc steps, # lambda points, vertices to change accuracy/runtime
steps1 = ceil(height/5) + 1; % discretization steps, each step is 5 nm
p1 = trirod(diameter, height, [12, 6, steps1]); 

steps2 = ceil(height10/5) + 1; % discretization steps, each step is 5 nm
p2 = trirod(diameter10, height10, [12, 6, steps2]); 

steps3 = ceil(height20/5) + 1; % discretization steps, each step is 5 nm
p3 = trirod(diameter20, height20, [12, 6, steps3]); 

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
exc = galerkin.planewave( [ 0, 0, 1 ], [ 1, 0, 0 ] );

%  light wavelength in vacuum
lambda = linspace( 560, 660, 30);
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

% Geometric cross section
radius = diameter/2;
cgeo = pi *(3/4 * radius^2 * (height - radius * 2/3))^(2/3);

% Extinction efficiencies
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

