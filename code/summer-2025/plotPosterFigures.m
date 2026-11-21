% triangle variables
triangleSus = [100.2237976
    106.3023001
    116.6130066
    123.198688
    139.7789121
    117.1448529
    110.0706235
    162.306428
    121.3547066
    188.0933553
    192.7459151
    128.961584
    176.8849313
    187.960425
    208.9086771
    182.6171875
    167.175293
    214.2789662
    207.902354
    212.5432357
    205.8105469
    219.7633222
    180.545041
    202.7587891
    176.5136719
    218.4918463
    201.8342495
    208.4882811
    225.6023034
    242.4281225
    218.57436
    210.7573748
    201.5380859
    250.2480879
    257.9862982
    248.5585749
    234.5581055
    268.4121132
    267.8833008
    251.2207031
    295.37201
    290.3137207
    289.8864746
    330.9020996
    347.3510742
    322.7092177
    363.1771859
    388.1531
    396.3321686
    426.8509746
    452.3765802];

triangleSub = [98.57150316
104.8743457
117.7938491
128.746596
146.7897072
125.8571446
121.4721501
134.6251667
129.5204312
145.391202
152.4533585
136.5737408
139.2587066
147.7134973
163.3380353
143.6767578
132.9345703
161.3295972
162.0386854
160.7967675
161.8652344
146.7008606
142.2622472
161.1328125
137.6953125
140.4355645
159.1783404
151.9937575
125.7569253
128.0095115
134.4218105
138.6939064
143.6157227
148.9382997
158.8676982
156.5385997
157.2265625
178.4962654
184.9975586
178.1616211
210.8001709
208.4960938
208.9233398
233.4899
240.9668
228.602998
246.9726056
259.7656
262.1087193
273.2658386
280.8747262];

triangleXa = [0.017413829
0.019802991
0.022951983
0.02729182
0.02729182
0.030141441
0.033655514
0.033655514
0.039850146
0.043889275
0.048839551
0.051758473
0.051758473
0.055048475
0.05685547
0.063065992
0.063065992
0.063065992
0.063065992
0.069104843
0.07079969
0.073817052
0.073817052
0.080695239
0.088986938
0.088986938
0.093806395
0.096957062
0.100827323
0.102100945
0.112004619
0.112004619
0.112004619
0.12257301
0.132580935
0.138963149
0.1511
0.162418821
0.183
0.183
0.2096
0.2175
0.2321
0.2452
0.2679
0.267949192
0.295365074
0.2954
0.316987298
0.334476874
0.34891526];

triangleR = 4 * sqrt(3) * triangleXa ./ (1 - 2 * triangleXa);

% disk variables

diskSus = [38.45214844
48.21777344
58.39324594
77.96272337
102.5390625
104.9657241
108.8879928
114.1747177
108.0322266
115.9667969
125.2441406
159.8775566
140.6622604
153.7230492
163.6110306
150.8789063
144.3321168
157.0717603
148.6821458
154.5681685
160.0341797
166.5097475
178.527832
156.7109078
170.2880859
170.6542969
158.0308035
177.6705623
198.5824034
245.7310259
233.5811615
283.1849
323.3516783
360.0766137];

diskSub = [25.02441406
33.56933594
42.81883836
65.62756598
77.51464844
78.23458314
79.93583381
83.71405005
76.29394531
82.39746094
83.49609375
82.84964114
76.71169639
80.23588657
88.00104707
79.83398438
76.51498318
85.32469422
80.3185299
85.29768437
91.55273438
101.7460465
113.8916016
92.9340899
107.2998047
108.0322266
98.76429588
113.7600541
129.1967511
150.5940557
144.0128326
172.3541
192.3229814
206.8837591];

diskXa = [0.055555556
0.071428571
0.1
0.1
0.125
0.136363636
0.147058824
0.15
0.166666667
0.166666667
0.166666667
0.1957
0.2
0.2059
0.2083
0.2143
0.2143
0.2222
0.2273
0.2333
0.25
0.2647
0.2778
0.2778
0.285714286
0.3
0.3
0.3125
0.3333
0.35
0.35
0.375
0.388888889
0.4];

diskR = 4 * diskXa ./ (1 - 2 * diskXa);

% plot suspension custom fit

modelTri = @(params, triangleR) (params(1) + params(2) * triangleR) ./ (params(3) + params(4) * triangleR);
modelDisk = @(params, diskR) (params(1) + params(2) * diskR) ./ (params(3) + params(4) * diskR);

guessTri = [100, 150, 0.3, 1];
guessDisk = [-10, 170, 0.7, 0.7];

options = optimoptions('lsqcurvefit', 'Display', 'off');

paramsTri = lsqcurvefit(modelTri, guessTri, triangleR, triangleSus, [], [], options);
paramsDisk = lsqcurvefit(modelDisk, guessDisk, diskR, diskSus, [], [], options);

triFit = modelTri(paramsTri, triangleR);
diskFit = modelDisk(paramsDisk, diskR);


% R^2 and r for the rational fits
SSres_tri  = sum( (triangleSus - triFit).^2 );
SStot_tri  = sum( (triangleSus - mean(triangleSus)).^2 );
R2_tri_rat = 1 - SSres_tri / SStot_tri;
r_tri_rat  = corr(triangleSus(:), triFit(:));  % Pearson r

SSres_disk  = sum( (diskSus - diskFit).^2 );
SStot_disk  = sum( (diskSus - mean(diskSus)).^2 );
R2_disk_rat = 1 - SSres_disk / SStot_disk;
r_disk_rat  = corr(diskSus(:), diskFit(:));    % Pearson r

% Strings to print on the figure with R and R²
eqTriStr  = sprintf('Triangle: (%.3f + %.3f R) / (%.3f + %.3f R)\nR = %.3f, R² = %.3f', ...
                    paramsTri(1), paramsTri(2), paramsTri(3), paramsTri(4), r_tri_rat, R2_tri_rat);

eqDiskStr = sprintf('Disk: (%.3f + %.3f R) / (%.3f + %.3f R)\nR = %.3f, R² = %.3f', ...
                    paramsDisk(1), paramsDisk(2), paramsDisk(3), paramsDisk(4), r_disk_rat, R2_disk_rat);

% Plot
figure();
plot(triangleR, triangleSus, 'ro', 'MarkerSize', 4, 'DisplayName', 'Triangle Data'); hold on
plot(triangleR, triFit,     'r-', 'LineWidth', 1.5, 'DisplayName', 'Triangle Fit');
plot(diskR,     diskSus,    'bo', 'MarkerSize', 4, 'DisplayName', 'Disk Data');
plot(diskR,     diskFit,    'b-', 'LineWidth', 1.5, 'DisplayName', 'Disk Fit');
legend('Location','best');
xlabel('Aspect Ratio R');
ylabel('Sensitivity S');
title('S vs R for Triangles and Disks in Suspension');

% Position text in the top-left corner
xLimits = xlim;
yLimits = ylim;
xTextPos = xLimits(1) + 0.02 * (xLimits(2) - xLimits(1));
yTextTri = yLimits(2) - 0.05 * (yLimits(2) - yLimits(1));
yTextDisk = yLimits(2) - 0.20 * (yLimits(2) - yLimits(1)); % Leave more room for multiline text

% Add text annotations
text(xTextPos, yTextTri, eqTriStr, 'Color', 'r', 'FontSize', 8, 'VerticalAlignment', 'top');
text(xTextPos, yTextDisk, eqDiskStr, 'Color', 'b', 'FontSize', 8, 'VerticalAlignment', 'top');

hold off



% plot substrate linear fit

pTri = polyfit(triangleXa, triangleSub, 1);
triFit = polyval(pTri, triangleXa);

pDisk = polyfit(diskXa, diskSub, 1);
diskFit = polyval(pDisk, diskXa);


% R^2 (and r) for linear fits
SSres_tri_lin  = sum( (triangleSub - triFit).^2 );
SStot_tri_lin  = sum( (triangleSub - mean(triangleSub)).^2 );
R2_tri_lin     = 1 - SSres_tri_lin / SStot_tri_lin;
r_tri_lin      = corr(triangleXa(:), triangleSub(:));   % Pearson r

SSres_disk_lin = sum( (diskSub - diskFit).^2 );
SStot_disk_lin = sum( (diskSub - mean(diskSub)).^2 );
R2_disk_lin    = 1 - SSres_disk_lin / SStot_disk_lin;
r_disk_lin     = corr(diskXa(:), diskSub(:));

% Equation strings: y = m x + b
if pTri(2) >= 0
    eqTriLin = sprintf('Triangle: y = %.3f x + %.3f', pTri(1), pTri(2));
else
    eqTriLin = sprintf('Triangle: y = %.3f x - %.3f', pTri(1), abs(pTri(2)));
end

if pDisk(2) >= 0
    eqDiskLin = sprintf('Disk: y = %.3f x + %.3f', pDisk(1), pDisk(2));
else
    eqDiskLin = sprintf('Disk: y = %.3f x - %.3f', pDisk(1), abs(pDisk(2)));
end


% Plot
figure();
plot(triangleXa, triangleSub, 'ro', 'MarkerSize', 4, 'DisplayName', 'Triangle Data'); hold on
plot(triangleXa, triFit,     'r-', 'LineWidth', 1.5, 'DisplayName', 'Triangle Linear Fit');
plot(diskXa,     diskSub,    'bo', 'MarkerSize', 4, 'DisplayName', 'Disk Data');
plot(diskXa,     diskFit,    'b-', 'LineWidth', 1.5, 'DisplayName', 'Disk Linear Fit');
legend('Location','best');
xlabel('Fractional Contact Area Xa');
ylabel('Sensitivity S');
title('S vs Xa for Triangles and Disks on Substrate');

% Annotate
text(0.02, 0.95, { ...
    sprintf('Triangle: y = %.3f x + %.3f', pTri(1), pTri(2)), ...
    sprintf('R^2 = %.4f  (r = %.4f)', R2_tri_lin, r_tri_lin) }, ...
    'Units','normalized','Color','r','VerticalAlignment','top','FontSize',10);

text(0.02, 0.78, { ...
    sprintf('Disk: y = %.3f x + %.3f', pDisk(1), pDisk(2)), ...
    sprintf('R^2 = %.4f  (r = %.4f)', R2_disk_lin, r_disk_lin) }, ...
    'Units','normalized','Color','b','VerticalAlignment','top','FontSize',10);

hold off
