% NT and ND on substrate lambda calculations

% Importing data from Excel file
% for Triangle, side x height
% for Disk, diameter x height
excelTriData = readcell('Simulation Data.xlsx', 'Sheet', 'Substrate % Reduction', 'Range', 'B2:B52');
excelDiskData = readcell('Simulation Data.xlsx', 'Sheet', 'Substrate % Reduction', 'Range', 'B60:B93');

numTri = numel(excelTriData);
triSide = zeros(1, numTri);
triHeight = zeros(1, numTri);

numDisk = numel(excelDiskData);
diskDiameter = zeros(1, numDisk);
diskHeight = zeros(1, numDisk);

for i = 1:numTri
    entry = excelTriData{i};
    entry = strrep(entry, 'x', ' x ');  % Add spaces around 'x' if missing
    parts = strsplit(entry, ' x '); % remove x from middle
    triSide(i) = str2double(parts{1}); % assign first value to side
    triHeight(i) = str2double(parts{2}); % assign second value to height
end

for j = 1:numDisk
    entry = excelDiskData{j};
    entry = strrep(entry, 'x', ' x ');  % Add spaces around 'x' if missing
    parts = strsplit(entry, ' x ');
    diskDiameter(i) = str2double(parts{1});
    diskHeight(i) = str2double(parts{2});
end

%% Using dimensions to calculate lambda
% suspension n = 1.33
% on substrate n = 1.52

nBulk = 1.33;
nSub = 1.52;

triLambdaSub = zeros(1, numTri);
triLambdaBulk = zeros(1, numTri);

diskLambdaSub = zeros(1, numDisk);
diskLambdaBulk = zeros(1, numDisk);

stepSize = 3;
refun = 3;
vertices = 20;

for k = 1:numTri
    triX0 = 500 + 4*k; % rough initial guess
    [triLambdaSub(k), triLambdaBulk(k)] = triangleLambdaMax(nSub, nBulk, triSide(k), triHeight(k), triX0, stepSize, refun);
end

for l = 1:numDisk
    diskX0 = 450 + 3*l; % rough initial guess
    [diskLambdaSub(l), diskLambdaBulk(l)] = diskLambdaMax(nSub, nBulk, diskDiameter(l), diskHeight(k), diskX0, stepSize, vertices);
end

