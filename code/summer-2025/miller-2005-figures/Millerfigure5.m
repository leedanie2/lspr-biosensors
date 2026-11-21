% Miller 2005 Figure 5 for Nanocylinders

% Dielectric properties
mat1a = Material( 1.33 ^ 2, 1 );
mat1b = Material( 1.41 ^ 2, 1 );
mat2 = Material( epstable( 'gold.dat' ), 1 );

%  material properties
mat_a = [ mat1a, mat2 ];
mat_b = [ mat1b, mat2 ];

% Nanocylinders (12 x 20, 30 x 20, 40 x 20, 50 x 20, 60 x 20)
poly = polygon(25, 'size', [20, 20]);

heights = [12, 30, 40, 50, 60];
rules = quadboundary.rules( 'quad3', triquad( 1 ) );
tau = cell(5, 2);
bem = cell(5, 2);

for i = 1:5
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
lambda = linspace( 500, 850, 100 );
k0 = 2 * pi ./ lambda;

%  far-field spectrometer
spec_a = galerkin.spectrum( mat_a, trisphere( 256 ), 'rules', rules );
spec_b = galerkin.spectrum( mat_b, trisphere( 256 ), 'rules', rules );

multiWaitbar( 'BEM solver', 0, 'Color', 'g', 'CanCancel', 'on' );

% Initialize before looping over wavenumbers
nk0 = numel(k0);
npol = size(exc.pol, 1);
cext12_a = zeros(nk0, npol );
cext12_b = zeros(nk0, npol );
cext30_a = zeros(nk0, npol );
cext30_b = zeros(nk0, npol );
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
  sol12_a = bem{1,1} \ exc( tau{1,1}, k0( i ) );
  sol12_b = bem{1,2} \ exc( tau{1,2}, k0( i ) );
  sol30_a = bem{2,1} \ exc( tau{2,1}, k0( i ) );
  sol30_b = bem{2,2} \ exc( tau{2,2}, k0( i ) );
  sol40_a = bem{3,1} \ exc( tau{3,1}, k0( i ) );
  sol40_b = bem{3,2} \ exc( tau{3,2}, k0( i ) );
  sol50_a = bem{4,1} \ exc( tau{4,1}, k0( i ) );
  sol50_b = bem{4,2} \ exc( tau{4,2}, k0( i ) );
  sol60_a = bem{5,1} \ exc( tau{5,1}, k0( i ) );
  sol60_b = bem{5,2} \ exc( tau{5,2}, k0( i ) );

  %  extinction cross sections
  cext12_a( i, : ) = extinction( exc, sol12_a );
  cext12_b( i, : ) = extinction( exc, sol12_b );
  cext30_a( i, : ) = extinction( exc, sol30_a );
  cext30_b( i, : ) = extinction( exc, sol30_b );
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

% geometric cross section
cgeo = pi * (3000 * heights).^(2/3);

% Preparing variables for plot
exteff12_a = cext12_a/cgeo(1);
exteff12_b = cext12_b/cgeo(1);
exteff30_a = cext30_a/cgeo(2);
exteff30_b = cext30_b/cgeo(2);
exteff40_a = cext40_a/cgeo(3);
exteff40_b = cext40_b/cgeo(3);
exteff50_a = cext50_a/cgeo(4);
exteff50_b = cext50_b/cgeo(4);
exteff60_a = cext60_a/cgeo(5);
exteff60_b = cext60_b/cgeo(5);
%% 

% Function for finding LSPR Wavelength
function lspr = findLSPR( wavelengths, exteff )
    [~, idx] = max(exteff);
    lspr = wavelengths(idx);
end

% Variables for plotting LSPR wavelength
ar = heights/20;
lspr12 = findLSPR(lambda, exteff12_a);
lspr30 = findLSPR(lambda, exteff30_a);
lspr40 = findLSPR(lambda, exteff40_a);
lspr50 = findLSPR(lambda, exteff50_a);
lspr60 = findLSPR(lambda, exteff60_a); 

%% 

% Function for finding sensitivity
function sens = findSens( wavelengths, n1, n2, exteff_1, exteff_2)
    lspr_1 = findLSPR(wavelengths, exteff_1);
    lspr_2 = findLSPR(wavelengths, exteff_2);
    sens = (lspr_2 - lspr_1)/(n2 - n1);
end

% Variables for plotting sensitivity
sens12 = findSens(lambda, 1.33, 1.41, exteff12_a, exteff12_b);
sens30 = findSens(lambda, 1.33, 1.41, exteff30_a, exteff30_b);
sens40 = findSens(lambda, 1.33, 1.41, exteff40_a, exteff40_b);
sens50 = findSens(lambda, 1.33, 1.41, exteff50_a, exteff50_b);
sens60 = findSens(lambda, 1.33, 1.41, exteff60_a, exteff60_b);

%% final plot for figure 5

figure(5);
plot(lspr12, sens12, 'o', 'MarkerEdgeColor','r', 'MarkerFaceColor','r'); hold on
plot(lspr30, sens30, 'o', 'MarkerEdgeColor','r', 'MarkerFaceColor','r'); hold on
plot(lspr40, sens40, 'o', 'MarkerEdgeColor','r', 'MarkerFaceColor','r'); hold on
plot(lspr50, sens50, 'o', 'MarkerEdgeColor','r', 'MarkerFaceColor','r'); hold on
plot(lspr60, sens60, 'o', 'MarkerEdgeColor','r', 'MarkerFaceColor','r'); hold off

xlabel( 'Lambda *' );
ylabel( 'RI Sensitivity' );

% Figure 5 Code for nanocylinders (double click to expand)