% Miller 2005 Figure 1 for Nanocylinders

% Dielectric properties
mat1a = Material( 1.33 ^ 2, 1 );
mat1b = Material( 1.41 ^ 2, 1 );
mat2 = Material( epstable( 'gold.dat' ), 1 );

%  material properties
mat_a = [ mat1a, mat2 ];
mat_b = [ mat1b, mat2 ];

% Nanocylinders (40 x 20, 50 x 20, 60 x 20)
poly = polygon(25, 'size', [20, 20]);

heights = [40, 50, 60];
rules = quadboundary.rules( 'quad3', triquad( 1 ) );
tau = cell(3, 2);
bem = cell(3, 2);

for i = 1:3
    edge = edgeprofile( heights(i), 9, 'mode', '11');
    p = tripolygon( poly, edge);
    p.verts(:,3) = p.verts(:,3) - min(p.verts(:,3)) + 1e-2;
    tau{i,1} = BoundaryEdge(mat_a, p, [1,2]);
    tau{i,2} = BoundaryEdge(mat_b, p, [1,2]);
    bem{i,1} = galerkin.bemsolver(tau{i,1}, 'rules', rules, 'waitbar', 1);
    bem{i,2} = galerkin.bemsolver(tau{i,2}, 'rules', rules, 'waitbar', 1);
end
%% 

% planewave excitation
exc = galerkin.planewave( [ 0, 0, 1 ], [ 1, 0, 0 ] );

%  light wavelength in vacuum
lambda = linspace( 600, 850, 100 );
k0 = 2 * pi ./ lambda;

%  far-field spectrometer
spec_a = galerkin.spectrum( mat_a, trisphere( 256 ), 'rules', rules );
spec_b = galerkin.spectrum( mat_b, trisphere( 256 ), 'rules', rules );

multiWaitbar( 'BEM solver', 0, 'Color', 'g', 'CanCancel', 'on' );

% Initialize before looping over wavenumbers
nk0 = numel(k0);
npol = size(exc.pol, 1);
cext40_a = zeros(nk0, npol );
cext40_b = zeros(nk0, npol );
cext50_a = zeros(nk0, npol );
cext50_b = zeros(nk0, npol );
cext60_a = zeros(nk0, npol );
cext60_b = zeros(nk0, npol );
%% 

%  loop over wavenumbers
for i = 1 : nk0
  %  solution of BEM equations
  sol40_a = bem{1,1} \ exc( tau{1,1}, k0( i ) );
  sol40_b = bem{1,2} \ exc( tau{1,2}, k0( i ) );

  sol50_a = bem{2,1} \ exc( tau{2,1}, k0( i ) );
  sol50_b = bem{2,2} \ exc( tau{2,2}, k0( i ) );

  sol60_a = bem{3,1} \ exc( tau{3,1}, k0( i ) );
  sol60_b = bem{3,2} \ exc( tau{3,2}, k0( i ) );

  %  extinction cross sections
  cext40_a( i, : ) = extinction( exc, sol40_a );
  cext40_b( i, : ) = extinction( exc, sol40_b );
  cext50_a( i, : ) = extinction( exc, sol50_a );
  cext50_b( i, : ) = extinction( exc, sol50_b );
  cext60_a( i, : ) = extinction( exc, sol60_a );
  cext60_b( i, : ) = extinction( exc, sol60_b );

  multiWaitbar( 'BEM solver', i / numel( k0 ) );
end
%% 

%  close waitbar
multiWaitbar( 'CloseAll' );

% Geometric cross section
cgeo = pi * (75 * heights).^(2/3);

% Preparing variables for plot
exteff40_a = cext40_a/cgeo(1);
exteff40_b = cext40_b/cgeo(1);
exteff50_a = cext50_a/cgeo(2);
exteff50_b = cext50_b/cgeo(2);
exteff60_a = cext60_a/cgeo(3);
exteff60_b = cext60_b/cgeo(3);

%%  final plot for Figure 1
figure(1);
plot(lambda, exteff40_a,'b'); hold on
plot(lambda, exteff40_b,'r'); hold on
plot(lambda, exteff50_a,'b--'); hold on
plot(lambda, exteff50_b,'r--'); hold on
plot(lambda, exteff60_a,'b:'); hold on
plot(lambda, exteff60_b,'r:'); hold off

xlabel( 'Wavelength, nm' );
ylabel( 'Extinction Efficiency' );
legend('40 x 20 n = 1.33', '40 x 20 n = 1.41', '50 x 20 n = 1.33', '50 x 20 n = 1.41', '60 x 20 n = 1.33', '60 x 20 n = 1.41');

% Figure 1 Code (double click to expand)