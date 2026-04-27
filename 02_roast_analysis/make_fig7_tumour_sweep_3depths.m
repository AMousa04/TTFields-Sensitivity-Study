% make_fig7_tumour_sweep_3depths.m
%
% Generates Figure 7 of the dissertation: percentage change in tumour ROI
% field metrics (Fan-like vs default conductivity) across five lateral
% positions at three axial depths.
%
% Input : Extended_TumourLocationSweep_3Depths.csv
% Output: Fig_TumourSweep_3Depths_PctChange.png
%
% The input CSV contains 15 ROI conditions (5 positions x 3 depths), each
% with:
%   names       - condition label, e.g. 'Xplus_Z83'
%   VoxelCount  - number of GM/WM voxels included in the sphere
%   A_median    - default-conductivity ROI median |E| (V/m)
%   B_median    - Fan-like ROI median |E| (V/m)
%   PctMedian   - percentage change in median (Fan-like vs default)
%   A_P95       - default ROI 95th-percentile |E| (V/m)
%   B_P95       - Fan-like ROI 95th-percentile |E| (V/m)
%   PctP95      - percentage change in P95 (Fan-like vs default)

clear; clc;

% ---- 1. Load data --------------------------------------------------------
csvFile = 'Extended_TumourLocationSweep_3Depths.csv';
T = readtable(csvFile);

% ---- 2. Parse position and depth from the condition name ----------------
parts      = split(T.names, '_');
T.Position = parts(:, 1);
T.Depth    = parts(:, 2);

% ---- 3. Set plotting order ----------------------------------------------
positions    = {'Base','Xplus','Xminus','Yplus','Yminus'};
posLabels    = {'Base','X+','X-','Y+','Y-'};
depths       = {'Z98','Z83','Z68'};
depthLabels  = {'z = 98 (superior)','z = 83 (reference)','z = 68 (inferior)'};

nPos   = numel(positions);
nDepth = numel(depths);

% ---- 4. Build matrices of percentage changes (rows = positions,
%        columns = depths) -------------------------------------------------
pctMedian = zeros(nPos, nDepth);
pctP95    = zeros(nPos, nDepth);

for i = 1:nPos
    for j = 1:nDepth
        idx = strcmp(T.Position, positions{i}) & strcmp(T.Depth, depths{j});
        pctMedian(i, j) = T.PctMedian(idx);
        pctP95(i, j)    = T.PctP95(idx);
    end
end

% ---- 5. Plot ------------------------------------------------------------
fig = figure('Color','w','Position',[100 100 1300 500]);

% Panel A: Median
subplot(1,2,1);
b1 = bar(pctMedian, 'grouped');
set(gca, 'XTickLabel', posLabels);
ylabel('% change in tumour ROI median |E| (Fan-like vs default)');
title('(A) Median');
yline(0, 'k-', 'LineWidth', 0.8, 'HandleVisibility','off');
legend(depthLabels, 'Location','southwest');
grid on;
set(gca, 'FontSize', 11);

% Panel B: P95
subplot(1,2,2);
b2 = bar(pctP95, 'grouped');
set(gca, 'XTickLabel', posLabels);
ylabel('% change in tumour ROI P95 |E| (Fan-like vs default)');
title('(B) 95th percentile');
yline(0, 'k-', 'LineWidth', 0.8, 'HandleVisibility','off');
legend(depthLabels, 'Location','northeast');
grid on;
set(gca, 'FontSize', 11);

% ---- 6. Save ------------------------------------------------------------
outFile = 'Fig_TumourSweep_3Depths_PctChange.png';
exportgraphics(fig, outFile, 'Resolution', 300);
fprintf('Saved: %s\n', outFile);
