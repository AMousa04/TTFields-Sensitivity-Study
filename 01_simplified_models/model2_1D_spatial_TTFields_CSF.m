%% Model 2 - 1D spatial TTFields with CSF (finite differences)
% Solve d/dx(sigma(x) dV/dx) = 0 along scalp -> skull -> CSF -> brain -> tumour.
% Continuous V(x) and E(x) along depth - useful sanity check before the 2D PDE
% model. CSF added explicitly because it dominates the field redistribution.

clear; clc; close all;

%% Geometry and grid
d_scalp  = 6e-3;
d_skull  = 7e-3;
d_csf    = 3e-3;
d_brain  = 20e-3;
d_tumour = 10e-3;

L_total = d_scalp + d_skull + d_csf + d_brain + d_tumour;

N  = 1000;
x  = linspace(0, L_total, N)';
dx = x(2) - x(1);

% layer boundary positions (end of each layer)
b1 = d_scalp;
b2 = b1 + d_skull;
b3 = b2 + d_csf;
b4 = b3 + d_brain;

%% Conductivity profile sigma(x)
% sigma values from McCann et al. 2019, CSF from Baumann et al. 1997
sigma = zeros(N,1);
for i = 1:N
    if      x(i) <= b1, sigma(i) = 0.41;   % scalp
    elseif  x(i) <= b2, sigma(i) = 0.02;   % skull
    elseif  x(i) <= b3, sigma(i) = 1.71;   % CSF
    elseif  x(i) <= b4, sigma(i) = 0.27;   % brain
    else                sigma(i) = 0.40;   % tumour
    end
end

%% Build finite-difference system
% Conservative form: sigma at the half-cell faces, harmonic-style averaging
A = sparse(N,N);
b = zeros(N,1);

for i = 2:N-1
    sL = (sigma(i) + sigma(i-1)) / 2;
    sR = (sigma(i) + sigma(i+1)) / 2;
    A(i,i-1) =  sL / dx^2;
    A(i,i)   = -(sL + sR) / dx^2;
    A(i,i+1) =  sR / dx^2;
end

%% Boundary conditions (Dirichlet at both ends)
V_left  =  1.0;
V_right = -1.0;

A(1,1) = 1;   b(1) = V_left;
A(N,N) = 1;   b(N) = V_right;

%% Solve
V = A \ b;
E = -gradient(V, dx);

%% Layer averages
idx_scalp  =  x <= b1;
idx_skull  = (x >  b1) & (x <= b2);
idx_csf    = (x >  b2) & (x <= b3);
idx_brain  = (x >  b3) & (x <= b4);
idx_tumour =  x >  b4;

E_scalp_avg  = mean(E(idx_scalp));
E_skull_avg  = mean(E(idx_skull));
E_csf_avg    = mean(E(idx_csf));
E_brain_avg  = mean(E(idx_brain));
E_tumour_avg = mean(E(idx_tumour));

%% Plots
figure;
plot(x*1000, sigma, 'LineWidth', 1.5);
xlabel('Depth from scalp (mm)'); ylabel('\sigma (S/m)');
title('Conductivity profile \sigma(x) (with CSF)');
grid on;

figure;
plot(x*1000, V, 'LineWidth', 1.5); hold on;
xline(b1*1000,'--k','Scalp/Skull');
xline(b2*1000,'--k','Skull/CSF');
xline(b3*1000,'--k','CSF/Brain');
xline(b4*1000,'--k','Brain/Tumour');
xlabel('Depth from scalp (mm)'); ylabel('V (V)');
title('Potential V(x)');
grid on;

figure;
plot(x*1000, E, 'LineWidth', 1.5); hold on;
xline(b1*1000,'--k'); xline(b2*1000,'--k');
xline(b3*1000,'--k'); xline(b4*1000,'--k');
xlabel('Depth from scalp (mm)'); ylabel('E (V/m)');
title('Electric field E(x)');
grid on;

figure;
E_layer_avg = [E_scalp_avg, E_skull_avg, E_csf_avg, E_brain_avg, E_tumour_avg];
bar(E_layer_avg);
set(gca,'XTickLabel',{'Scalp','Skull','CSF','Brain','Tumour'},'FontSize',12);
ylabel('Average E (V/m)');
title('Layer-averaged E-field (1D spatial model)');
grid on;

%% Summary
fprintf('Model 2 - 1D spatial TTFields with CSF\n\n');
fprintf('Mean |E| per layer (V/m):\n');
fprintf('  Scalp  : %.3f\n', E_scalp_avg);
fprintf('  Skull  : %.3f\n', E_skull_avg);
fprintf('  CSF    : %.3f\n', E_csf_avg);
fprintf('  Brain  : %.3f\n', E_brain_avg);
fprintf('  Tumour : %.3f\n', E_tumour_avg);