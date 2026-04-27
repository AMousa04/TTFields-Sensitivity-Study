%% Model 1 - 1D series-circuit TTFields (scalp -> skull -> brain -> tumour)
% Sanity-check forward model. Same J in every layer, E_layer = J / sigma.
% Used to build intuition before moving to 2D / ROAST.

clear; clc; close all;

%% TTFields settings
f       = 200e3;           % 200 kHz
V0      = 1.0;             % 1 V across the stack for now
nCycles = 5;
Ttot    = nCycles / f;
Fs      = 20 * f;          % 20 samples / cycle
t       = 0 : 1/Fs : Ttot;

%% Layer thicknesses (m) and conductivities (S/m)
% sigma values from McCann et al. 2019 review at ~200 kHz, brain lumped
% between GM 0.47 and WM 0.22.
d.scalp  = 6e-3;
d.skull  = 7e-3;
d.brain  = 20e-3;
d.tumour = 10e-3;

sigma.scalp  = 0.41;
sigma.skull  = 0.02;
sigma.brain  = 0.27;
sigma.tumour = 0.40;       % tumour slightly more conductive than normal brain

layerNames = {'Scalp','Skull','Brain','Tumour'};
d_vals     = [d.scalp,     d.skull,     d.brain,     d.tumour];
sigma_vals = [sigma.scalp, sigma.skull, sigma.brain, sigma.tumour];

%% Solve series circuit
% V0 = J * sum(d_i / sigma_i)  =>  J = V0 / sum(d_i / sigma_i)
denom = sum(d_vals ./ sigma_vals);
J     = V0 / denom;            % A/m^2 (uniform along the stack)
E_vals = J ./ sigma_vals;      % V/m in each layer

E_scalp  = E_vals(1);
E_tumour = E_vals(4);

%% Time-domain waveforms at scalp and tumour
E_scalp_t  = E_scalp  * sin(2*pi*f*t);
E_tumour_t = E_tumour * sin(2*pi*f*t);

%% Plots
figure;
bar(E_vals);
set(gca,'XTickLabel',layerNames,'FontSize',12);
ylabel('E-field magnitude (V/m)');
title('TTFields E-field per layer (1D series model)');
grid on;

figure;
plot(t*1e6, E_scalp_t,  'LineWidth', 1.5); hold on;
plot(t*1e6, E_tumour_t, 'LineWidth', 1.5);
xlabel('Time (\mus)');
ylabel('E-field (V/m)');
legend('Scalp','Tumour','Location','best');
title('TTFields waveform: scalp vs tumour');
grid on;

%% Summary
fprintf('Model 1 - 1D series conduction\n');
fprintf('f = %.1f kHz, V0 = %.2f V\n\n', f/1e3, V0);
for i = 1:numel(layerNames)
    fprintf('%-7s : d = %4.1f mm, sigma = %.3f S/m, E = %.3f V/m\n', ...
        layerNames{i}, d_vals(i)*1e3, sigma_vals(i), E_vals(i));
end
fprintf('\nE_tumour          = %.3f V/m\n',  E_tumour);
fprintf('E_tumour/E_scalp  = %.3f\n',         E_tumour/E_scalp);