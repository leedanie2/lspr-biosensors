% Parameters
vertices = 20;
diameter = 40;
analyteThickness = 10;
diameterOuter = diameter + 2 * analyteThickness;
height = 12;
stepSize = 5;
steps = ceil(height / stepSize) + 1;
edge = edgeprofile(height, steps, 'mode', '11');

% Outer polygon: positive direction (normal outward)
polyOuter = polygon(vertices, 'size', [diameterOuter, diameterOuter], 'dir', +1);
% Inner polygon (hole): negative direction (normal inward)
polyInner = polygon(vertices, 'size', [diameter, diameter], 'dir', -1);
% Combine polygons into array to represent polygon with hole
polyRing = [polyOuter, polyInner];
% Extrude into 3D ring shape
pRing = tripolygon(polyRing, edge);

% % Plot to verify
% plot(pRing);
% axis equal;
% title('3D Nanodisk Ring Geometry with Hole');

% Make cap
polyCap = polygon(vertices, 'size', [diameterOuter, diameterOuter], 'dir', +1);
heightCap = analyteThickness;
stepsCap = ceil(heightCap / stepSize) + 1;
edgeCap = edgeprofile(heightCap, stepsCap, 'mode', '11');
pCap = tripolygon(polyCap, edgeCap);

% Shift cap upward by height of ring
pCap = shift(pCap, [0, 0, height]);

% % Plot full analyte layer
% figure;
% plot(pRing); hold on
% plot(pCap); hold off
% axis equal;
% view(3);
% title('Complete Analyte Layer with Ring and Cap');

% Metal nanoparticle
poly = polygon(vertices, 'size', [diameter, diameter]);
p = tripolygon(poly, edge);

% Plot all three components
figure;
plot(pRing, 'FaceColor', [0 0.7 1], 'FaceAlpha', 0.4);   % cyan-ish ring
hold on;
plot(pCap,  'FaceColor', [1 1 0],   'FaceAlpha', 0.6);   % yellow cap
plot(p,'FaceColor', [1 0 0],   'FaceAlpha', 0.7);   % red metal

axis equal;
title('Two-Layer System on Substrate: Metal + Analyte Cap + Ring');

