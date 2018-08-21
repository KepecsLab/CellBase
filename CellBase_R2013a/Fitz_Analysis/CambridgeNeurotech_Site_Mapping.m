%% E2 Probe
% see documentation for ADPT-HS36 adapter and cambridge neurotech probes...
neuralynx = (1:32)';
neuralynx_zeroBased = neuralynx - 1;
cambridge = ([11 9 7 5 3 1 2 4 6 8 10 12 13 14 15 16 ...
    17 18 19 20 21 23 22 25 24 27 26 29 28 31 30 32])';

pins = ([1 2 3 4 5 15 6 16 7 8 9 10 17 18 19 20 ... % pin IDs from ADPT-HS36
    30 29 28 27 40 39 38 37 26 36 25 35 34 33 32 31])';

E2_trodes   = ([16 12 8 6 2 3 5 9 15 10 14 13 4 1 7 11 ...
    17 21 22 25 27 29 28 30 18 23 19 20 24 26 31 32])';



probe_mapping = table(pins, cambridge, neuralynx, neuralynx_zeroBased);

% example- access neuralynx channel numbers corresponding to E2 probe
% trodes
table_sorted = sortrows(probe_mapping, 'cambridge'); % order table by cambridge indices
table_sorted(E2_trodes, :)  % display rows comprising (in order) TT1 - TT8

%% P2 Probe

P2_trodes   = ([16 3 12 5 8 9 6 11 2 7 4 1 14 13 15 10 ...
    17 29 21 28 22 30 25 32 27 31 24 26 19 20 18 23])';



probe_mapping_P2 = table(pins, cambridge, neuralynx, neuralynx_zeroBased);

% example- access neuralynx channel numbers corresponding to E2 probe
% trodes
table_sorted_P2 = sortrows(probe_mapping, 'cambridge'); % order table by cambridge indices
table_sorted_P2(P2_trodes, :)  % display rows comprising (in order) TT1 - TT8

%% make P2 geom.csv file and define trode ncs file number to ms4 trode site number mapping

% put cambridge neurotech site IDs into 4 columns where columns 1:2
% correspond to shank A, and coumns 3:4 correspond to shank B
P2_trodes_ms4 = [16 12 8 6 2 4 14 15 3 5 9 11 7 1 13 10 ...
    17 21 22 25 27 24 19 18 29 28 30 32 31 26 20 23];
P2_trodes_ms4 = P2_trodes_ms4(:);
P2_trodes_ms4 = reshape(P2_trodes_ms4, 8, 4);

% make 8 x 2 x 2 matrix (8 sites by 2 trodes by x,y) to hold x,y coordinates of each site on a shank
P2_geo = zeros(8,2,2);
% x coordinates
P2_geo(:,1,1) = 0; P2_geo(:,2,1) = 22.5; % 22.5um spacing between columns
% y coordinates
P2_geo(:,:,2) = repmat((0:25:(25*(8-1)))', 1, 2); P2_geo(:,2,2) = P2_geo(:,2,2) + 12.5; % shift second row by 12.5um in y


%% make table for ms4 loading engine and geom.csv file creation
trodeIx = [ones(16,1); (ones(16,1) + 1)];
trodeSiteIx = (1:32)'; 
x = P2_geo(:,:,1); x = x(:); x = repmat(x, 2, 1);
y = P2_geo(:,:,2); y = y(:); y = repmat(y, 2, 1);
P2_ms4_geo = table(trodeIx, trodeSiteIx, table_sorted_P2.neuralynx(P2_trodes_ms4(:)),...
    x, y, 'VariableNames', {'trode', 'site', 'ncs', 'x', 'y'});

%% make geom.csv file (nsites x (x,y))
P2_ms4_geo = sortrows(P2_ms4_geo, 'site');

% write the csv file
toCsv = zeros(16,2);
toCsv(:,1) = P2_ms4_geo.x(1:16);
toCsv(:,2) = P2_ms4_geo.y(1:16);

[fname, pname] = uiputfile('path', 'Choose geometry file base directory');
pname = [pname 'P2' filesep];
ensureDirectory(pname);
csvwrite(fullfile(pname, 'geo.csv'), toCsv);
disp(['*** saving geo.csv in ' pname ' ***']);

%% make trode creation file mapping .mat file (nSites x nTrodes)
P2_ms4_geo = sortrows(P2_ms4_geo, {'trode', 'site'});
toMat = zeros(16,2);
toMat(:) = P2_ms4_geo.ncs;
save(fullfile(pname, 'ncs2ms4_map.mat'), 'toMat');