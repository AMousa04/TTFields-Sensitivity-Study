function make_bundle_pack(runID, projectFolder, labelText)
% MAKE_BUNDLE_PACK  Bundle one ROAST run into a tidy reproducibility folder.
%   make_bundle_pack(runID, projectFolder, labelText)
%
%   Copies the key raw outputs from rawDir into a dated run folder, runs the
%   GM/WM/(GM+WM) metric extraction, and writes axial+coronal |E| slice PNGs
%   plus a brain |E| histogram. Used to keep the two final dissertation runs
%   (174510 and 182002) side by side without manually shuffling files.
%
%   Example:
%   make_bundle_pack('20260211T182002','01_Model_Replication_FanInspired','Fan cond + FT7->T8')

if nargin < 1, error('runID required'); end
if nargin < 2, projectFolder = '01_Model_Replication_FanInspired'; end
if nargin < 3, labelText = ''; end

base   = 'C:\ROAST\roast\example';
rawDir = fullfile(base, projectFolder, 'run_outputs');

%% Folder structure for this run bundle
destRoot  = fullfile(base, projectFolder, 'runs', ['run_' runID]);
inputsDir = fullfile(destRoot, 'inputs');
outDir    = fullfile(destRoot, 'outputs_raw');
anaDir    = fullfile(destRoot, 'analysis');
figDir    = fullfile(destRoot, 'figures');

if ~exist(inputsDir,'dir'), mkdir(inputsDir); end
if ~exist(outDir,   'dir'), mkdir(outDir);    end
if ~exist(anaDir,   'dir'), mkdir(anaDir);    end
if ~exist(figDir,   'dir'), mkdir(figDir);    end

%% Copy the raw run outputs
% copyIf is just a "copy if it exists" shortcut so the script doesn't fall
% over when an optional file isn't present (e.g. .pos files only exist for
% some runs).
copyIf = @(f,dst) (exist(f,'file')==2) && copyfile(f,dst);

emag = fullfile(rawDir, ['MNI152_T1_1mm_' runID '_emag.nii']);   % |E|
vnii = fullfile(rawDir, ['MNI152_T1_1mm_' runID '_v.nii']);      % V

copyIf(emag, outDir);
copyIf(vnii, outDir);

% Optional ROAST artefacts - copied if present, ignored if not
optMat   = fullfile(rawDir, ['MNI152_T1_1mm_' runID '_roastOpt.mat']);
resMat   = fullfile(rawDir, ['MNI152_T1_1mm_' runID '_roastResult.mat']);
usedElec = fullfile(rawDir, ['MNI152_T1_1mm_' runID '_usedElec.mat']);
msh      = fullfile(rawDir, ['MNI152_T1_1mm_' runID '.msh']);
pro      = fullfile(rawDir, ['MNI152_T1_1mm_' runID '.pro']);
epos     = fullfile(rawDir, ['MNI152_T1_1mm_' runID '_e.pos']);
vpos     = fullfile(rawDir, ['MNI152_T1_1mm_' runID '_v.pos']);
readyMsh = fullfile(rawDir, ['MNI152_T1_1mm_' runID '_ready.msh']);
maskElec = fullfile(rawDir, ['MNI152_T1_1mm_' runID '_mask_elec.nii']);
maskGel  = fullfile(rawDir, ['MNI152_T1_1mm_' runID '_mask_gel.nii']);
enii     = fullfile(rawDir, ['MNI152_T1_1mm_' runID '_e.nii']);

copyIf(optMat,   outDir);
copyIf(resMat,   outDir);
copyIf(usedElec, outDir);
copyIf(msh,      outDir);
copyIf(pro,      outDir);
copyIf(epos,     outDir);
copyIf(vpos,     outDir);
copyIf(readyMsh, outDir);
copyIf(maskElec, outDir);
copyIf(maskGel,  outDir);
copyIf(enii,     outDir);

%% Reproducibility inputs (segmentation mask, log, RAS-aligned MRI)
roastLog = fullfile(rawDir, 'MNI152_T1_1mm_roastLog');
maskFile = fullfile(rawDir, 'MNI152_T1_1mm_ras_T1orT2_SPM_masks.nii');
rasMRI   = fullfile(rawDir, 'MNI152_T1_1mm_ras.nii');

copyIf(roastLog, inputsDir);
copyIf(maskFile, inputsDir);
copyIf(rasMRI,   inputsDir);

%% Pack folder for derived outputs (metrics CSV + 3 figures)
packDir = fullfile(destRoot, ['pack_' runID]);
if ~exist(packDir,'dir'), mkdir(packDir); end

if exist(emag,'file') ~= 2
    error('Cannot find emag: %s', emag);
end
if exist(maskFile,'file') ~= 2
    error('Cannot find mask file: %s', maskFile);
end

E = double(niftiread(emag));
M = double(niftiread(maskFile));

% SPM12 tissue labels via ROAST: 1=GM, 2=WM, 3=CSF, 4=bone, 5=skin, 6=air
GM    = (M==1);
WM    = (M==2);
BRAIN = GM | WM;

%% Metrics: Nvox, mean, median, P10, P90, P95, P99, max for GM, WM, GM+WM
getMetrics = @(x) [numel(x), mean(x), median(x), ...
                   prctile(x,10), prctile(x,90), prctile(x,95), prctile(x,99), ...
                   max(x)];

% Drop non-finite and zero voxels (zeros sit outside the head)
Egm = E(GM);    Egm = Egm(isfinite(Egm) & Egm>0);
Ewm = E(WM);    Ewm = Ewm(isfinite(Ewm) & Ewm>0);
Ebr = E(BRAIN); Ebr = Ebr(isfinite(Ebr) & Ebr>0);

T = array2table([getMetrics(Egm); getMetrics(Ewm); getMetrics(Ebr)], ...
    'VariableNames', {'Nvox','Mean','Median','P10','P90','P95','P99','Max'}, ...
    'RowNames',      {'GM','WM','Brain_GMplusWM'});

disp(T);
writetable(T, fullfile(packDir,'Metrics_GM_WM_Brain.csv'), 'WriteRowNames', true);
copyfile(fullfile(packDir,'Metrics_GM_WM_Brain.csv'), anaDir);

%% Slice picks - axial and coronal slice with the most brain voxels
% (i.e. the slice that gives the cleanest visual)
brainCount_ax  = squeeze(sum(sum(BRAIN, 2), 3));
[~, sAx]       = max(brainCount_ax);

brainCount_cor = squeeze(sum(sum(BRAIN, 1), 3));
[~, sCor]      = max(brainCount_cor);

%% Axial slice plot
imgAx = squeeze(E(sAx,:,:));
imgAx(~squeeze(BRAIN(sAx,:,:))) = NaN;       % mask out non-brain so colorbar stays useful

figure('Color','w');
imagesc(imgAx); axis image off; colorbar;
set(gca,'YDir','normal');
title(sprintf('|E| axial (run %s)', runID), 'Interpreter','none');
exportgraphics(gcf, fullfile(packDir,'Axial_Efield.png'), 'Resolution', 300);

%% Coronal slice plot
imgCor = squeeze(E(:,sCor,:));
imgCor(~squeeze(BRAIN(:,sCor,:))) = NaN;

figure('Color','w');
imagesc(imgCor'); axis image off; colorbar;
set(gca,'YDir','normal');
title(sprintf('|E| coronal (run %s)', runID), 'Interpreter','none');
exportgraphics(gcf, fullfile(packDir,'Coronal_Efield.png'), 'Resolution', 300);

%% Brain |E| histogram (right-skewed - motivates reporting median + P95)
figure('Color','w');
histogram(Ebr, 100);
xlabel('|E| in brain (V/m)');
ylabel('Voxel count');
title(sprintf('Brain |E| distribution (run %s)', runID), 'Interpreter','none');
exportgraphics(gcf, fullfile(packDir,'BrainHistogram.png'), 'Resolution', 300);

copyfile(fullfile(packDir,'Axial_Efield.png'),    figDir);
copyfile(fullfile(packDir,'Coronal_Efield.png'),  figDir);
copyfile(fullfile(packDir,'BrainHistogram.png'),  figDir);

%% Manifest - small text file recording what was bundled
manifest = fullfile(destRoot, 'manifest.txt');
fid = fopen(manifest, 'w');
fprintf(fid, 'Run ID: %s\n',     runID);
fprintf(fid, 'Label: %s\n\n',    labelText);
fprintf(fid, 'MATLAB: %s\n',     version);
fprintf(fid, 'ROAST:  %s\n\n',   which('roast.m'));
fprintf(fid, 'RawDir: %s\n\n',   rawDir);
fprintf(fid, 'Core outputs:\n%s\n%s\n\n', emag, vnii);
fclose(fid);

fprintf('\nBundle: %s\nPack:   %s\n', destRoot, packDir);
end