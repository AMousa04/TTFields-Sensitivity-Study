%% Model 3b - Tumour position sweep in 2D
% Repeat Model 3 with tumour shifted along x. Quick check that ROI
% position genuinely affects |E| under this geometry, before committing
% to the 3D ROAST sweep where each run is much more expensive.

clear; clc; close all;

R_head   = 0.09;
R_tumour = 0.01;

% x positions to test (m): left of centre, centre, right of centre
tumour_positions = [-0.03, 0.00, 0.03];
tumour_cy = 0.0;

sigma_brain  = 0.27;
sigma_tumour = 0.40;

fprintf('Tumour position sweep (x in m)\n');
fprintf('------------------------------------------------\n');
fprintf('  x_centre   E_brain_avg   E_tumour_avg   ratio\n');

for k = 1:numel(tumour_positions)

    tumour_cx = tumour_positions(k);

    % build geometry for this position
    model = createpde();
    gd = [ 1      1;
           0      tumour_cx;
           0      tumour_cy;
           R_head R_tumour ];
    ns = char('H1','T1'); ns = ns';
    sf = 'H1+T1';

    [g,~] = decsg(gd, sf, ns);
    geometryFromEdges(model, g);
    generateMesh(model, 'Hmax', 0.005);

    % Face 1 = brain, Face 2 = tumour
    specifyCoefficients(model,'m',0,'d',0,'c',sigma_brain, 'a',0,'f',0,'Face',1);
    specifyCoefficients(model,'m',0,'d',0,'c',sigma_tumour,'a',0,'f',0,'Face',2);

    % same electrode convention as Model 3
    applyBoundaryCondition(model,'dirichlet','Edge',[1 4],'u', 1.0);
    applyBoundaryCondition(model,'dirichlet','Edge',[2 3],'u',-1.0);

    % solve
    result = solvepde(model);
    [gradx, grady] = evaluateGradient(result);
    E_mag = sqrt((-gradx).^2 + (-grady).^2);

    % element-wise mean |E|, then split by face
    brain_elems  = findElements(model.Mesh,'region','Face',1);
    tumour_elems = findElements(model.Mesh,'region','Face',2);
    E_elem       = mean(E_mag(model.Mesh.Elements), 1);

    E_brain_avg  = mean(E_elem(brain_elems));
    E_tumour_avg = mean(E_elem(tumour_elems));
    ratio        = E_tumour_avg / E_brain_avg;

    fprintf('%9.3f   %10.3f   %12.3f   %5.3f\n', ...
            tumour_cx, E_brain_avg, E_tumour_avg, ratio);

    % only plot the last case (saves a load of figure windows)
    if k == numel(tumour_positions)
        figure;
        pdeplot(model,'XYData',E_mag,'Contour','on');
        axis equal; colorbar;
        xlabel('x (m)'); ylabel('y (m)');
        title(sprintf('|E| for tumour at x = %.3f m', tumour_cx));
    end
end