function tumour_roi_sweep_3D(emag_default_file, emag_fan_file, mask_file, outDir)
% TUMOUR_ROI_SWEEP_3D  Spherical-ROI tumour metric sweep across positions/depths.
%
%   Loads two ROAST |E| volumes (default and Fan-like conductivity) plus the
%   SPM12 segmentation mask, builds a spherical ROI of radius 20 voxels at
%   each of 5 lateral positions and 3 axial depths, intersects each sphere
%   with the GM/WM mask, and extracts Nvox/Mean/Median/P10/P90/P95/Max from
%   |E| inside the ROI.
%
%   Writes three CSVs to outDir:
%     - TumourLocationSweep_<runID>_default.csv
%     - TumourLocationSweep_<runID>_fan.csv
%     - Extended_TumourLocationSweep_3Depths.csv  (paired comparison)
%
%   ROI definition is sphere ∩ (GM | WM). Full sphere alone did not
%   reproduce the original z=83 reference spreadsheet - see logbook
%   Section 8 for the debugging story.
%
%   Example call:
%     tumour_roi_sweep_3D( ...
%       'MNI152_T1_1mm_20260211T174510_emag.nii', ...
%       'MNI152_T1_1mm_20260211T182002_emag.nii', ...
%       'MNI152_T1_1mm_ras_T1orT2_SPM_masks.nii', ...
%       'tables_out');

if nargin < 4, outDir = pwd; end
if ~exist(outDir,'dir'), mkdir(outDir); end

%% Load volumes
assert(isfile(emag_default_file), 'Default emag not found: %s', emag_default_file);
assert(isfile(emag_fan_file),     'Fan emag not found: %s',     emag_fan_file);
assert(isfile(mask_file),         'Mask file not found: %s',    mask_file);

E_default = double(niftiread(emag_default_file));
E_fan     = double(niftiread(emag_fan_file));
M         = double(niftiread(mask_file));

assert(isequal(size(E_default), size(M)), 'default emag and mask differ');
assert(isequal(size(E_fan),     size(M)), 'fan emag and mask differ');

% SPM12 labels: 1=GM, 2=WM, 3=CSF, 4=bone, 5=skin, 6=air
brainMask = (M==1) | (M==2);

%% ROI grid - 5 positions x 3 depths = 15 conditions
% Reference centre at [93,104,83]; lateral offsets ±25 vox; depths z=98,83,68
r        = 20;                              % radius in voxels (= mm; isotropic 1 mm)
base_xy  = [93 104];
offsets  = [  0   0;                        % Base
             25   0;                        % X+
            -25   0;                        % X-
              0  25;                        % Y+
              0 -25];                       % Y-
posNames = {'Base','Xplus','Xminus','Yplus','Yminus'};
depths   = [98 83 68];

% Pre-compute coordinate grids once (faster than re-meshing inside the loop)
[X,Y,Z] = ndgrid(1:size(M,1), 1:size(M,2), 1:size(M,3));

%% Sweep
nROI = numel(posNames) * numel(depths);
rows_default = cell(nROI, 1);
rows_fan     = cell(nROI, 1);
rows_pair    = cell(nROI, 1);

idx = 0;
for di = 1:numel(depths)
    cz = depths(di);
    for oi = 1:size(offsets,1)
        idx = idx + 1;
        cx = base_xy(1) + offsets(oi,1);
        cy = base_xy(2) + offsets(oi,2);

        % sphere ∩ brain (GM or WM only)
        sphere    = (X-cx).^2 + (Y-cy).^2 + (Z-cz).^2 <= r^2;
        tumourROI = sphere & brainMask;
        Nvox      = nnz(tumourROI);

        % extract |E| values inside the ROI for both runs
        Ed = E_default(tumourROI);
        Ef = E_fan(tumourROI);

        % per-run metrics
        rows_default{idx} = makeRow(posNames{oi}, cx, cy, cz, Nvox, Ed);
        rows_fan{idx}     = makeRow(posNames{oi}, cx, cy, cz, Nvox, Ef);

        % paired comparison row (matches Extended_TumourLocationSweep_3Depths.csv)
        med_d = median(Ed,'omitnan');
        med_f = median(Ef,'omitnan');
        p95_d = prctile(Ed, 95);
        p95_f = prctile(Ef, 95);
        rows_pair{idx} = { ...
            sprintf('%s_Z%d', posNames{oi}, cz), ...
            Nvox, ...
            med_d, med_f, 100*(med_f - med_d)/med_d, ...
            p95_d, p95_f, 100*(p95_f - p95_d)/p95_d };
    end
end

%% Write per-run CSVs
% Layout matches TumourLocationSweep_20260211T*.csv:
%   Location,CenterX,CenterY,CenterZ,Nvox,Mean,Median,P10,P90,P95,Max
T_default = cell2table(vertcat(rows_default{:}), 'VariableNames', ...
    {'Location','CenterX','CenterY','CenterZ','Nvox','Mean','Median','P10','P90','P95','Max'});
T_fan = cell2table(vertcat(rows_fan{:}), 'VariableNames', ...
    {'Location','CenterX','CenterY','CenterZ','Nvox','Mean','Median','P10','P90','P95','Max'});

writetable(T_default, fullfile(outDir, 'TumourLocationSweep_default.csv'));
writetable(T_fan,     fullfile(outDir, 'TumourLocationSweep_fan.csv'));

%% Write paired-comparison CSV
% Layout matches Extended_TumourLocationSweep_3Depths.csv:
%   names,VoxelCount,A_median,B_median,PctMedian,A_P95,B_P95,PctP95
T_pair = cell2table(vertcat(rows_pair{:}), 'VariableNames', ...
    {'names','VoxelCount','A_median','B_median','PctMedian','A_P95','B_P95','PctP95'});

writetable(T_pair, fullfile(outDir, 'Extended_TumourLocationSweep_3Depths.csv'));

fprintf('\nWrote 3 CSVs to %s\n', outDir);
disp(T_pair);
end


%% ---- local helper: build one CSV row from |E| values inside a ROI ----
function row = makeRow(name, cx, cy, cz, Nvox, Evals)
    row = { name, cx, cy, cz, Nvox, ...
            mean(Evals,'omitnan'), ...
            median(Evals,'omitnan'), ...
            prctile(Evals,10), ...
            prctile(Evals,90), ...
            prctile(Evals,95), ...
            max(Evals) };
end