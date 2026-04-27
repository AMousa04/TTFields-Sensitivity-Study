%% Model 3 (layered) - 2D TTFields in a layered head with embedded tumour
% Layers (outside in): scalp -> skull -> CSF -> brain -> tumour.
% Solve ∇·(sigma(x,y)·∇V) = 0 with Dirichlet ±1 V on opposite scalp arcs.
% Builds on Model 3: same 2D head, but with proper layered conductivity
% structure rather than uniform brain. CSF in particular dominates how
% current is redistributed before reaching the tumour.

clear; clc; close all;

%% Geometry
R_head = 0.09;        % 90 mm

% layer thicknesses (m)
d_scalp = 6e-3;
d_skull = 7e-3;
d_csf   = 3e-3;

% inner-boundary radii from centre, working inwards
R_skull_outer = R_head        - d_scalp;     % scalp/skull
R_csf_outer   = R_skull_outer - d_skull;     % skull/CSF
R_brain       = R_csf_outer   - d_csf;       % CSF/brain

% tumour: 10 mm radius, offset 30 mm right
R_tumour  = 0.01;
tumour_cx = 0.03;
tumour_cy = 0.00;

%% Conductivities (S/m)
sigma_scalp  = 0.41;
sigma_skull  = 0.02;
sigma_csf    = 1.71;
sigma_brain  = 0.27;
sigma_tumour = 0.40;

%% Build geometry (head + tumour as union, layers handled in coefficient fn)
model = createpde();

gd = [ 1      1;
       0      tumour_cx;
       0      tumour_cy;
       R_head R_tumour ];
ns = char('H1','T1'); ns = ns';
sf = 'H1+T1';

[dl,~] = decsg(gd, sf, ns);
geometryFromEdges(model, dl);

figure;
pdegplot(model,'EdgeLabels','on');
axis equal;
title('Layered head geometry (outer head + tumour)');

%% Mesh
generateMesh(model, 'Hmax', 0.004);    % ~4 mm
figure;
pdemesh(model);
axis equal;
title('FE mesh (layered head)');

%% Spatially varying sigma(x,y) - via local function below
cFunc = @(location,state) sigmaLayered( ...
    location, ...
    R_head, R_skull_outer, R_csf_outer, R_brain, ...
    tumour_cx, tumour_cy, R_tumour, ...
    sigma_scalp, sigma_skull, sigma_csf, sigma_brain, sigma_tumour);

specifyCoefficients(model, 'm',0,'d',0,'c',cFunc,'a',0,'f',0);

%% Boundary conditions (same convention as Model 3)
edges_pos = [1 4];
edges_neg = [2 3];

applyBoundaryCondition(model,'dirichlet','Edge',edges_pos,'u', 1.0);
applyBoundaryCondition(model,'dirichlet','Edge',edges_neg,'u',-1.0);

%% Solve
result = solvepde(model);
V = result.NodalSolution;

[gradx, grady] = evaluateGradient(result);
Ex = -gradx;
Ey = -grady;
E_mag = sqrt(Ex.^2 + Ey.^2);

%% Plots
figure;
pdeplot(model,'XYData',V,'Contour','on','ColorMap','parula');
axis equal; xlabel('x (m)'); ylabel('y (m)');
title('V(x,y) in layered head');
colorbar;

figure;
pdeplot(model,'XYData',E_mag,'Contour','on');
axis equal; xlabel('x (m)'); ylabel('y (m)');
title('|E(x,y)| in layered head');
colorbar;

%% Layer-averaged |E|
% Compute element centroids and average |E| per element. Then assign each
% element to a tissue based on its radius from head centre (or distance
% to tumour centre for the tumour mask).
p        = model.Mesh.Nodes;
elements = model.Mesh.Elements(1:3,:);     % triangles only - first 3 rows
nElem    = size(elements,2);

xc     = zeros(1,nElem);
yc     = zeros(1,nElem);
E_elem = zeros(1,nElem);

for e = 1:nElem
    n        = elements(:,e);
    xc(e)    = mean(p(1,n));
    yc(e)    = mean(p(2,n));
    E_elem(e) = mean(E_mag(n));
end

r_elem = sqrt(xc.^2 + yc.^2);
r_tum  = sqrt((xc - tumour_cx).^2 + (yc - tumour_cy).^2);

idx_scalp  = (r_elem > R_skull_outer) & (r_elem <= R_head);
idx_skull  = (r_elem > R_csf_outer)   & (r_elem <= R_skull_outer);
idx_csf    = (r_elem > R_brain)       & (r_elem <= R_csf_outer);
idx_brain  =  r_elem <= R_brain;       % includes tumour - excluded below
idx_tumour =  r_tum  <= R_tumour;
idx_brain_only = idx_brain & ~idx_tumour;

E_scalp_avg  = mean(E_elem(idx_scalp));
E_skull_avg  = mean(E_elem(idx_skull));
E_csf_avg    = mean(E_elem(idx_csf));
E_brain_avg  = mean(E_elem(idx_brain_only));
E_tumour_avg = mean(E_elem(idx_tumour));

%% Bar plot
figure;
bar([E_scalp_avg, E_skull_avg, E_csf_avg, E_brain_avg, E_tumour_avg]);
set(gca,'XTickLabel',{'Scalp','Skull','CSF','Brain','Tumour'},'FontSize',12);
ylabel('Average |E| (V/m)');
title('Layer-averaged |E| (2D layered model)');
grid on;

%% Summary
fprintf('\nModel 3 (layered) - 2D layered head + tumour\n\n');
fprintf('  scalp  : %.3f V/m\n', E_scalp_avg);
fprintf('  skull  : %.3f V/m\n', E_skull_avg);
fprintf('  CSF    : %.3f V/m\n', E_csf_avg);
fprintf('  brain  : %.3f V/m\n', E_brain_avg);
fprintf('  tumour : %.3f V/m\n', E_tumour_avg);
fprintf('  tumour/brain ratio : %.3f\n', E_tumour_avg/E_brain_avg);


%% ---- local function: sigma(x,y) ----
function c = sigmaLayered(location, ...
    R_head, R_skull_outer, R_csf_outer, R_brain, ...
    tumour_cx, tumour_cy, R_tumour, ...
    sigma_scalp, sigma_skull, sigma_csf, sigma_brain, sigma_tumour)

x = location.x;
y = location.y;

r   = sqrt(x.^2 + y.^2);
r_t = sqrt((x - tumour_cx).^2 + (y - tumour_cy).^2);

% start with brain everywhere, overwrite by layer working outwards,
% then overwrite tumour region last
c = sigma_brain * ones(size(r));
c(r > R_brain       & r <= R_csf_outer)   = sigma_csf;
c(r > R_csf_outer   & r <= R_skull_outer) = sigma_skull;
c(r > R_skull_outer & r <= R_head)        = sigma_scalp;
c(r_t <= R_tumour) = sigma_tumour;
end