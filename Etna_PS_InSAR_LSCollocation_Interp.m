clc
clear
close all

addpath('SAR_Data');
addpath('Matlab_Functions');

% --------------------- LEAST SQUARE COLLOCATION ------------------------ %

%% 1. SAR data 2D - Cleaned input for collocation
% Import raw SAR data
data = readmatrix('PS_Etna_SE_new.csv');
easting = data(2:end,2);
northing = data(2:end,3);
proj_from = projcrs(3035);
[lat, lon] = projinv(proj_from, easting, northing);

mean_vel = data(2:end,4); % mean_velocity (unit: mm/year)
mean_vel_std = data(2:end,5); % mean_velocity_standard_deviation (mm/year)
displ2D = data(2:end, 6:end);
day_rel = data(1,6:end);

% Convert latitude and longitude into UTM x/y
[x, y] = deg2utm(lat, lon);
[~, idx_sort] = sort(x);
x_sort = x(idx_sort);
y_sort = y(idx_sort);

% === Load processed spline outputs ===
load('Etna_spline_output_cleaned.mat') 
% includes: x_sort3, y_sort3, z_in, y_est_fin, best_res_epc, displ_fin, res_poly, Cyy_est

% Use cleaned and sorted data (after outlier removal)
x_sort = x_sort3;
y_sort = y_sort3;
displ_epc = z_in;               % mean velocity after outlier removal
poly = y_est_fin;               % polynomial trend re-estimated
residuals_poly = res_poly;      % displacement - trend
poly_spl = displ_fin;           % spline model = poly + residual
poly_std = sqrt(diag(Cyy_est)); % standard deviation of polynomial surface

% ===== Plot 1: Zero-Mean Residuals Distribution =====
figure
plot3((x_sort-mean(x_sort))/1000, (y_sort-mean(y_sort))/1000, res_poly, 'x')
grid on
xlabel('X (km)');
ylabel('Y (km)');
zlabel('Residual Mean Velocity (mm/year)');
title('Zero-Mean Residuals Distribution');

% ===== Plot 2: 3D Scatter - Original, Polynomial, Spline+Poly =====
figure
scatter3(x_sort, y_sort, displ_epc, 10, 'filled')
hold on
scatter3(x_sort, y_sort, poly, 10, 'filled')
scatter3(x_sort, y_sort, poly_spl, 10, 'filled')
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Mean Velocity (mm/year)');
title('3D View: Observed vs Polynomial vs Spline Surface');
legend('Observed (mean velocity)', 'Polynomial Trend', 'Spline + Polynomial');

% ===== Plot 3: 2D comparison - Polynomial and Residuals =====
figure

subplot(1,2,1)
scatter(x_sort - mean(x_sort), y_sort - mean(y_sort), 10, poly, 'filled');
xlabel('X (centered) (m)'); ylabel('Y (centered) (m)');
title('2D Polynomial Trend Surface');
colorbar

subplot(1,2,2)
scatter(x_sort - mean(x_sort), y_sort - mean(y_sort), 10, poly_spl - poly, 'filled');
xlabel('X (centered) (m)'); ylabel('Y (centered) (m)');
title('Spline Residuals (Spline - Poly)');
colorbar

% Residuals used for collocation
residuals_poly = displ_epc - poly;

%% 2. Empirical covariance
% Compute empirical covariance for 2D spatial residuals using f2DCovEmpEst.
% The sample distance is used as the representative distance for each bin.
[sigma, eCov, sigmaGrid, eCovF, Cecf, ~] = f2DCovEmpEst(residuals_poly, x_sort, y_sort, 100, 'cartesian', 2);

% ===== Plot 4: Empirical Covariance Curve =====
figure
plot(sigmaGrid, eCovF, '.-', 'LineWidth', 1.3, 'MarkerSize', 12);
xlabel('Distance (m)');
ylabel('Empirical Covariance (mm^2/year^2)');
title('Empirical Covariance vs Distance');
xlim([0 16000]);
grid on

% ===== Plot 5: Covariance Weights =====
figure
plot(sigmaGrid, sqrt(1 ./ Cecf), 'LineWidth', 1.5);
xlabel('Distance (m)');
ylabel('Weight (1/sqrt(Var))');
title('Weighting Function for Covariance Fitting');
grid on

%% 3. Covariance modelling
idx_zero = find(eCovF < 0, 1);  % first point less than 0
if isempty(idx_zero)
    last_bin = length(eCovF);   % no negative value, keep all
else
    last_bin = idx_zero - 1;    % use the last bin before the negative value
end

% Exponential model
% a parameter
a = ((log(eCovF(2))-log(eCovF(3)))/(sigmaGrid(3) - sigmaGrid(2)));

% k parameter
Q = diag(Cecf(2:17));
yo = log(eCovF(2:17)) - log(4.2) - log(cos(pi*sigmaGrid(2:17)/3000));
A = -a * sigmaGrid(2:17);
% k could be estimated by GLS; here we choose it manually:
k = 0.9;

fcove = @(tau) 4.2 * exp(-a * k * tau) .* cos(pi * 1/3000 * tau);

% Plot model vs empirical
figure
plot(sigmaGrid, eCovF, '.-', 'LineWidth', 1.3, 'MarkerSize', 12);
title('Covariance Modelling');
xlim([0 16000]);
hold on;
plot(sigmaGrid, fcove(sigmaGrid), 'LineWidth', 2);

% Check the "nugget" amplitude (compare with observation error level)
sv_exp = sqrt(eCovF(1) - fcove(0));
fprintf('\nNugget from Exponential modelling %.4f\n', sv_exp);

% A-priori (0.3 mm)
s2v_ap = 0.3^2;
fprintf('A-priori noise amplitude (std) %.4f\n', sqrt(s2v_ap));

% Selected covariance model
fcov = fcove;
s2v = sv_exp^2;

%% 4. Collocation - Filtering
% Compute covariance matrices and apply collocation for 2D spatial residuals.

% Estimation locations (same as observations)
x_estim = x_sort;
y_estim = y_sort;

% Estimation–observation covariance (Cv_vo)
[x1, x2] = meshgrid(x_estim, x_sort);
[y1, y2] = meshgrid(y_estim, y_sort);
dist_vo = sqrt((x1 - x2).^2 + (y1 - y2).^2);
Cv_vo = fcov(dist_vo);
figure; imagesc(dist_vo); title('Distance Matrix (est–obs)'); colorbar

% Observation–observation covariance (Cvo_vo_filt)
[x3, x4] = meshgrid(x_sort, x_sort);
[y3, y4] = meshgrid(y_sort, y_sort);
dist_vo_vo = sqrt((x3 - x4).^2 + (y3 - y4).^2);
Cvo_vo_filt = fcov(dist_vo_vo);

% Collocation estimate at observation points (filtering)
v_est =  Cv_vo * inv(Cvo_vo_filt + s2v .* eye(length(x_sort))) * residuals_poly;
Cy0y0 = Cvo_vo_filt + s2v .* eye(length(x_sort));
displ_coll  = poly + v_est; 

% Plot residuals vs collocation estimate
figure;
scatter3(x_sort, y_sort, residuals_poly, 10, 'b', 'filled'); % original residuals
hold on;
scatter3(x_sort, y_sort, v_est, 10, 'r', 'filled');   % collocation estimate
title('residuals vs collocation estimate');
legend('Residuals (obs - poly)','Collocation estimate (v_{est})');

%% 5. Collocation - Prediction
% Predict residuals on a regular grid and reconstruct the field.

% Regular grid for estimation
xGrid = min(x):100:max(x);
yGrid = min(y):100:max(y);
[xGridMesh, yGridMesh] = meshgrid(xGrid, yGrid);
x_est_grid = xGridMesh(:);  % flatten
y_est_grid = yGridMesh(:);

% Estimation–observation covariance for grid
[x_obs, x_est_p] = meshgrid(x_sort, x_est_grid);
[y_obs, y_est_p] = meshgrid(y_sort, y_est_grid);
dist_vo_grid = sqrt((x_obs - x_est_p).^2 + (y_obs - y_est_p).^2);
Cv_vo_grid = fcov(dist_vo_grid);

% Collocation prediction on grid
v_est_grid = Cv_vo_grid * inv(Cvo_vo_filt + s2v .* eye(length(x_sort))) * residuals_poly;
displ_coll_grid = reshape(v_est_grid, size(xGridMesh));

%% 6. Error of the prediction
% Posterior standard deviation of the collocated residuals
W = (Cvo_vo_filt + s2v * eye(length(x_sort))) \ eye(length(x_sort));
v_var = fcov(0) - diag(Cv_vo * W * Cv_vo');
v_std = sqrt(v_var);

%% Residual Comparison
% 1) Residuals: collocation vs spline
res_coll = v_est;
res_spl  = poly_spl - poly;             % spline residuals

rms_fun  = @(x) sqrt(mean(x.^2));       % RMS

% Basic stats
stats = @(x) [mean(x), std(x), rms_fun(x), min(x), max(x)];
S_coll = stats(res_coll);
S_spl  = stats(res_spl);

% Pairwise comparison
diff_rs   = res_coll - res_spl;
corr_rs   = corr(res_coll(:), res_spl(:), 'rows','complete');
rms_diff  = rms_fun(diff_rs);
mae_diff  = mean(abs(diff_rs));
mxdiff    = max(abs(diff_rs));

fprintf('\n--- Residuals (Collocation) ---\n');
fprintf('mean=%.3f, std=%.3f, RMS=%.3f, min=%.3f, max=%.3f (mm/yr)\n', S_coll);
fprintf('--- Residuals (Spline) ---\n');
fprintf('mean=%.3f, std=%.3f, RMS=%.3f, min=%.3f, max=%.3f (mm/yr)\n', S_spl);

fprintf('\n[Colloc vs Spline residuals]\n');
fprintf('Corr=%.3f, RMS diff=%.3f, MAE=%.3f, Max|diff|=%.3f (mm/yr)\n', ...
        corr_rs, rms_diff, mae_diff, mxdiff);

% 2) Uncertainty summary (v_std)
p = prctile(v_std,[5 50 95]);
fprintf('\n--- Posterior Std v_std (collocation) ---\n');
fprintf('mean=%.3f, median=%.3f, p05=%.3f, p95=%.3f, min=%.3f, max=%.3f (mm/yr)\n', ...
        mean(v_std), p(2), p(1), p(3), min(v_std), max(v_std));

% 3) Final signal comparison: (poly + v_est) vs (poly + spline) = poly_spl
final_coll = displ_coll;    % poly + v_est
final_spl  = poly_spl;      % poly + (poly_spl - poly)

corr_final   = corr(final_coll(:), final_spl(:), 'rows','complete');
rms_final    = rms_fun(final_coll - final_spl);
maxabs_final = max(abs(final_coll - final_spl));

fprintf('\n[Final signal: Colloc (poly+v) vs Spline (poly+spline)]\n');
fprintf('Corr=%.3f, RMS diff=%.3f, Max|diff|=%.3f (mm/yr)\n', ...
        corr_final, rms_final, maxabs_final);

% 4) Normalized error (using dynamic range of spline residuals)
dyn_spl = max(res_spl) - min(res_spl);
NRMSE   = rms_diff / dyn_spl;
fprintf('\nNRMSE (residuals, vs spline range)=%.3f\n', NRMSE);

% 5) Combine polynomial & collocation uncertainty in quadrature
% (more appropriate than simple summation)
if exist('poly_std','var')
    tot_std_q = sqrt(v_std.^2 + poly_std.^2);
    pt = prctile(tot_std_q,[5 50 95]);
    fprintf('\n--- Total Std (quadrature) ---\n');
    fprintf('mean=%.3f, median=%.3f, p05=%.3f, p95=%.3f (mm/yr)\n', ...
            mean(tot_std_q), pt(2), pt(1), pt(3));
end

%% 7. Final signal
% Reconstruct final signal as polynomial trend + collocated residuals.

v_est_grid_interp = griddata(x_sort, y_sort, v_est, xGridMesh, yGridMesh, 'cubic');
v_est_grid_interp_rs = v_est_grid_interp(:);

% Plot - Filtering results and their std
figure
subplot(1,2,1)
scatter(x_sort-mean(x_sort), y_sort-mean(y_sort), 20, v_est, 'filled'); 
title('Collocation Estimated Residuals');
xlabel('X (m, zero-mean)'); ylabel('Y (m, zero-mean)');
colorbar; colormap turbo; grid on;

subplot(1,2,2)
scatter(x_sort-mean(x_sort), y_sort-mean(y_sort), 20, v_std, 'filled'); 
title('Std. Dev. of Collocation (v_{std})');
xlabel('X (m, zero-mean)'); ylabel('Y (m, zero-mean)');
colorbar; colormap turbo; grid on;

% Plot - Collocation vs Splines
figure
subplot(1,2,1)
scatter(x_sort-mean(x_sort), y_sort-mean(y_sort), 20, v_est, 'filled'); 
title('Collocation Residual Estimate');
xlabel('X (m, zero-mean)'); ylabel('Y (m, zero-mean)');
colorbar; colormap turbo;  grid on;

subplot(1,2,2)
scatter(x_sort-mean(x_sort), y_sort-mean(y_sort), 20, poly_spl - poly, 'filled'); 
title('Spline Residuals');
xlabel('X (m, zero-mean)'); ylabel('Y (m, zero-mean)');
colorbar; colormap turbo; grid on;

% Plot - Final estimate (filtering + polynomial)
figure
subplot(1,2,1)
scatter(x_sort-mean(x_sort), y_sort-mean(y_sort), 20, displ_coll, 'filled'); 
title('Final Estimated Signal (Poly + v)');
xlabel('X (m, zero-mean)'); ylabel('Y (m, zero-mean)');
colorbar; colormap turbo; grid on;

subplot(1,2,2)
scatter(x_sort-mean(x_sort), y_sort-mean(y_sort), 20, v_std + poly_std, 'filled'); 
title('Total Std ( sqrt( v_{std}^2 + poly_{std}^2 ) )');
xlabel('X (m, zero-mean)'); ylabel('Y (m, zero-mean)');
colorbar; colormap turbo; grid on;

if ~exist('xGridMesh','var') || ~exist('yGridMesh','var')
    xGrid = min(x_sort):100:max(x_sort);
    yGrid = min(y_sort):100:max(y_sort);
    [xGridMesh, yGridMesh] = meshgrid(xGrid, yGrid);
end

figure
% Scatter of the “observation velocity” (cleaned mean velocity)
scatter3(x_sort, y_sort, displ_epc, 10, 'filled'); 
hold on
% Interpolate (poly + v_est) onto the mesh and render as a surface
surf(xGridMesh, yGridMesh, ...
     griddata(x_sort, y_sort, displ_coll, xGridMesh, yGridMesh, 'cubic'), ...
     'EdgeColor','none','FaceAlpha',0.7);

xlabel('X (m)');
ylabel('Y (m)');
zlabel('Mean Velocity (mm/year)');
title('Final Collocation Surface (Poly + v)');
colormap(turbo); colorbar; grid on; view(3);

%% 8. Least Squares Collocation
% Design matrix for a 2nd-degree polynomial
A = [ones(size(x_sort)), x_sort, y_sort, x_sort.^2, x_sort.*y_sort, y_sort.^2];

% LS solution
x_est = (A' * inv(Cy0y0) * A) \ (A' * inv(Cy0y0) * displ_epc);
y_est = A * x_est;
res1 = displ_epc - y_est;
s02_est = res1' * inv(Cy0y0) * res1 / (size(A,1) - size(A,2));
Cxx_est = s02_est * inv(A' * inv(Cy0y0) * A);

% --- 8.1 Significance test on parameters
alpha = 0.05;
t_x_est = x_est ./ sqrt(diag(Cxx_est));
t_lim = tinv(1 - alpha/2, size(A,1) - size(A,2));

fprintf('\n\nt-student test on the significance of the estimated parameters (it. 0): \n')
for i = 1:length(x_est)
    fprintf('Param. %i, t_obs: %.4f (t_lim: %.4f)\n', i, abs(t_x_est(i)), t_lim);
end

poly_LSC = A * x_est;
res_LSC = displ_epc - poly_LSC;
% Difference between trends and RMS
diff_LSC = poly_LSC - poly;
diff_LSC_rms = rms(diff_LSC);

% Plot
figure

% --- LSC Trend Estimate ---
subplot(1,3,1)
scatter(x_sort - mean(x_sort), y_sort - mean(y_sort), 15, poly_LSC, 'filled')
title('LSC Trend Estimate')
xlabel('X (m, zero-mean)'); ylabel('Y (m, zero-mean)');
colorbar; colormap turbo; grid on;

% --- LS Trend Estimate ---
subplot(1,3,2)
scatter(x_sort - mean(x_sort), y_sort - mean(y_sort), 15, poly, 'filled')
title('LS Trend Estimate')
xlabel('X (m, zero-mean)'); ylabel('Y (m, zero-mean)');
colorbar; colormap turbo; grid on;

% --- Difference ---
subplot(1,3,3)
scatter(x_sort - mean(x_sort), y_sort - mean(y_sort), 15, diff_LSC, 'filled')
title('Trend Difference (LSC - LS)')
xlabel('X (m, zero-mean)'); ylabel('Y (m, zero-mean)');
colorbar; colormap turbo; grid on;

%% Iteration 1
% Covariance model from LSC residuals
[sigma2, eCov2, sigmaGrid2, eCovF2, Cecf2, ~] = f2DCovEmpEst(res_LSC, x_sort, y_sort, 100, 'cartesian', 0);

% Plot empirical covariance
figure
plot(sigmaGrid2, eCovF2, '.-', 'LineWidth', 1.3, 'MarkerSize', 12);
xlim([0 16000]);
title('Empirical Covariance');
xlabel('\tau (distance in meters)');
ylabel('Covariance');
grid on;

% Exponential model fit with cosine
a = ((log(eCovF2(2))-log(eCovF2(3)))/(sigmaGrid2(3) - sigmaGrid2(2)));
k = 0.8;
fcovLSC = @(tau) 5.7 * exp(-a * k * tau) .* cos(pi * 1/5500 * tau);
s2v_LSC = eCovF2(1) - fcovLSC(0);

% Plot model vs empirical
figure;
plot(sigmaGrid2, eCovF2, 'b.-', 'LineWidth', 1.3, 'MarkerSize', 10); hold on;
plot(sigmaGrid2, fcovLSC(sigmaGrid2), 'r-', 'LineWidth', 2);
xlim([0 16000]);
ylim([-4, max(eCovF2)+1]);
xlabel('\tau (distance in meters)');
ylabel('Covariance');
title('Covariance Model Fitting: Exponential Decay + Cosine');
legend('Empirical Covariance', 'Exponential + Cosine Model', 'Location', 'northeast');
grid on;

% Check the "nugget" amplitude (compare with observation error)
sv_exp2 = sqrt(eCovF2(1) - fcovLSC(0));
fprintf('\nNugget from Exponential modelling %.4f\n', sv_exp2);

% A-priori (0.3 mm)
s2v_ap_2 = 0.3^2;
fprintf('A-priori noise amplitude (std) %.4f\n', sqrt(s2v_ap));

% Selected covariance model for iteration 1
fcov2 = fcovLSC;
s2v = sv_exp^2;

%%
n = length(x_sort);
x0 = (x_sort - mean(x_sort))/1000;   % km
y0 = (y_sort - mean(y_sort))/1000;   % km
A_LSC2 = [ones(size(x0)), x0, y0, x0.^2, x0.*y0, y0.^2];

% === 8.3 Remodel covariance by res_LSC ===
% Ensure non-negativity for noise variance
s2v_LSC2 = max(0, eCovF(1) - fcovLSC(0));

% New observation–observation covariance and GLS weight matrix
Cvo_vo_LSC = fcovLSC(dist_vo_vo);
W_LSC      = inv(Cvo_vo_LSC + s2v_LSC2*eye(n));   % GLS weight matrix

% Iterate once: re-estimate trend under the new weight matrix
x_est_LSC2 = (A_LSC2' * W_LSC * A_LSC2) \ (A_LSC2' * W_LSC * displ_epc);
poly_LSC2  = A_LSC2 * x_est_LSC2;
res_LSC2   = displ_epc - poly_LSC2;

% GLS residual variance and parameter covariance
s02_est_LSC2 = res_LSC2' * W_LSC * res_LSC2 / (size(A_LSC2,1) - size(A_LSC2,2));
Cxx_est_LSC2 = s02_est_LSC2 * inv(A_LSC2' * W_LSC * A_LSC2);

% --- Significance test on parameters (LSC, it. 1)
alpha = 0.05;
t_x_est_LSC2 = x_est_LSC2 ./ sqrt(diag(Cxx_est_LSC2));
t_lim = tinv(1 - alpha/2, size(A_LSC2,1) - size(A_LSC2,2));

fprintf('\n\nt-student test on the significance of the estimated parameters (LSC, it. 1): \n')
for i = 1:length(x_est_LSC2)
    fprintf('Param. %i, t_obs: %.4f (t_lim: %.4f)\n', i, abs(t_x_est_LSC2(i)), t_lim);
end

% Comparison (RMS)
diff_LSC2 = poly_LSC2 - poly;
diff_LSC2_rms = rms(diff_LSC2);

rms_diff_poly_LSC = rms(poly_LSC2 - poly_LSC);
fprintf('\n[RMS] Trend difference (Iter.1 - Iter.0): %.4f mm/yr\n', rms_diff_poly_LSC);

%%
% 8.4 - New collocation with LSC
Cv_vo_LSC  = fcovLSC(dist_vo);                           % est–obs covariance
v_est_LSC  = Cv_vo_LSC * W_LSC * res_LSC2;               % filtering (same-point estimation)
displ_coll_LSC1 = poly_LSC2 + v_est_LSC;                 % final field (iteration 1)

% A posteriori variance (filtered)
v_var_LSC  = fcovLSC(0) - diag(Cv_vo_LSC * W_LSC * Cv_vo_LSC');
v_std_LSC  = sqrt(v_var_LSC);

% Compare LSC polynomial iterations
diff_poly_LSC = poly_LSC2 - poly_LSC;
rms_diff_poly_LSC = rms(diff_poly_LSC);
fprintf('\n[RMS] Trend difference between LSC_iter_1 and LSC_iter_0: %.4f mm/year\n', rms_diff_poly_LSC);

% --- 8.5 Plot
figure('Color','w');

% LSC trend (iteration 0)
subplot(1,3,1)
scatter(x_sort - mean(x_sort), y_sort - mean(y_sort), 15, poly_LSC, 'filled');
title('LSC Trend Estimate (Iter. 0)');
xlabel('X (m, zero-mean)'); ylabel('Y (m, zero-mean)');
colorbar; colormap turbo; grid on;

% LSC trend (iteration 1)
subplot(1,3,2)
scatter(x_sort - mean(x_sort), y_sort - mean(y_sort), 15, poly_LSC2, 'filled');
title('LSC Trend Estimate (Iter. 1)');
xlabel('X (m, zero-mean)'); ylabel('Y (m, zero-mean)');
colorbar; colormap turbo; grid on;

% Difference
subplot(1,3,3)
scatter(x_sort - mean(x_sort), y_sort - mean(y_sort), 15, diff_poly_LSC, 'filled');
title('Trend Difference (Iter.1 - Iter.0)');
xlabel('X (m, zero-mean)'); ylabel('Y (m, zero-mean)');
colorbar; colormap turbo; grid on;

% Comparison of estimates
figure('Color', 'w');
subplot(1,2,1)
scatter(x_sort - mean(x_sort), y_sort - mean(y_sort), 20, v_est_LSC, 'filled');
title('LSC Collocation Estimate v_{est}^{LSC}');
xlabel('X (m, zero-mean)');
ylabel('Y (m, zero-mean)');
colorbar; colormap turbo; grid on;

subplot(1,2,2)
scatter(x_sort - mean(x_sort), y_sort - mean(y_sort), 20, v_std_LSC, 'filled');
title('LSC Estimation Std. Dev. \sigma_{v}');
xlabel('X (m, zero-mean)');
ylabel('Y (m, zero-mean)');
colorbar; colormap turbo; grid on;

% Comparison with non-LSC collocation
figure('Color','w');
subplot(1,2,1)
scatter(x_sort - mean(x_sort), y_sort - mean(y_sort), 20, v_est, 'filled');
title('Standard Collocation Estimate v_{est}');
xlabel('X (m, zero-mean)');
ylabel('Y (m, zero-mean)');
colorbar; colormap turbo; grid on;

subplot(1,2,2)
scatter(x_sort - mean(x_sort), y_sort - mean(y_sort), 20, v_est_LSC, 'filled');
title('LSC Collocation Estimate v_{est}^{LSC}');
xlabel('X (m, zero-mean)');
ylabel('Y (m, zero-mean)');
colorbar; colormap turbo; grid on;

%% 9) Final LSC Fit
% Final LSC: polynomial (Iter.1) + LSC residuals
displ_LSC_final = poly_LSC2 + v_est_LSC;      % [mm/yr] at observation points

% Polynomial covariance on each point: sigma_poly_LSC2(i) = sqrt( a_i^T * Cxx * a_i )
poly_var_LSC2 = sum((A_LSC2 * Cxx_est_LSC2) .* A_LSC2, 2);
poly_std_LSC2 = sqrt(max(poly_var_LSC2, 0));  % numerically safe

% Total uncertainty
tot_std_LSC = sqrt(poly_std_LSC2.^2 + v_std_LSC.^2);

% --- Comparison with other methods
% 1) vs OLS+Collocation final field displ_coll
corr_LSC_vs_Coll = corr(displ_LSC_final(:), displ_coll(:), 'rows', 'complete');
rms_LSC_vs_Coll  = sqrt(mean((displ_LSC_final - displ_coll).^2));
maxabs_LSC_vs_Coll = max(abs(displ_LSC_final - displ_coll));

% 2) vs Spline final field poly_spl
corr_LSC_vs_Spl = corr(displ_LSC_final(:), poly_spl(:), 'rows', 'complete');
rms_LSC_vs_Spl  = sqrt(mean((displ_LSC_final - poly_spl).^2));
maxabs_LSC_vs_Spl = max(abs(displ_LSC_final - poly_spl));

% 3) Residuals of LSC final field vs observed mean velocity
residuals_LSC_final = displ_epc - displ_LSC_final;
rms_res_LSC_final   = sqrt(mean(residuals_LSC_final.^2));
mean_res_LSC_final  = mean(residuals_LSC_final);
std_res_LSC_final   = std(residuals_LSC_final);

fprintf('\n=== LSC Final Fit (at observations) ===\n');
fprintf('Residuals to observed mean-vel: mean=%.3f, std=%.3f, RMS=%.3f (mm/yr)\n', ...
        mean_res_LSC_final, std_res_LSC_final, rms_res_LSC_final);
fprintf('[LSC vs OLS+Collocation] Corr=%.3f, RMS=%.3f, Max|diff|=%.3f (mm/yr)\n', ...
        corr_LSC_vs_Coll, rms_LSC_vs_Coll, maxabs_LSC_vs_Coll);
fprintf('[LSC vs Spline]         Corr=%.3f, RMS=%.3f, Max|diff|=%.3f (mm/yr)\n', ...
        corr_LSC_vs_Spl,  rms_LSC_vs_Spl,  maxabs_LSC_vs_Spl);

% --- Build final LSC surface on a regular grid
if ~exist('xGridMesh','var') || ~exist('yGridMesh','var')
    xGrid = min(x_sort):100:max(x_sort);
    yGrid = min(y_sort):100:max(y_sort);
    [xGridMesh, yGridMesh] = meshgrid(xGrid, yGrid);
end

xg0 = (xGridMesh(:) - mean(x_sort)) / 1000;  
yg0 = (yGridMesh(:) - mean(y_sort)) / 1000;
A_grid = [ones(size(xg0)), xg0, yg0, xg0.^2, xg0.*yg0, yg0.^2];
poly_grid_LSC = reshape(A_grid * x_est_LSC2, size(xGridMesh));

% 2) LSC residual prediction on the grid
[x_obs_g, x_est_g] = meshgrid(x_sort, xGridMesh(:));
[y_obs_g, y_est_g] = meshgrid(y_sort, yGridMesh(:));
dist_vo_grid_LSC = sqrt((x_obs_g - x_est_g).^2 + (y_obs_g - y_est_g).^2);
Cv_vo_grid_LSC   = fcovLSC(dist_vo_grid_LSC);

v_est_grid_LSC   = Cv_vo_grid_LSC * W_LSC * res_LSC2;       % grid residuals
v_est_grid_LSC   = reshape(v_est_grid_LSC, size(xGridMesh));

% 3) Final grid field
displ_LSC_grid   = poly_grid_LSC + v_est_grid_LSC;

% --- Visualization (final field and uncertainty at observation points)
figure('Color','w');
subplot(1,3,1)
scatter(x_sort - mean(x_sort), y_sort - mean(y_sort), 20, displ_LSC_final, 'filled');
title('LSC Final Signal (obs points)');
xlabel('X (m, zero-mean)'); ylabel('Y (m, zero-mean)'); grid on; colorbar; colormap turbo;

subplot(1,3,2)
scatter(x_sort - mean(x_sort), y_sort - mean(y_sort), 20, tot_std_LSC, 'filled');
title('LSC Total Std (poly \oplus v)'); 
xlabel('X (m, zero-mean)'); ylabel('Y (m, zero-mean)'); grid on; colorbar; colormap turbo;

subplot(1,3,3)
scatter(x_sort - mean(x_sort), y_sort - mean(y_sort), 20, residuals_LSC_final, 'filled');
title('Residuals: observed - LSC final');
xlabel('X (m, zero-mean)'); ylabel('Y (m, zero-mean)'); grid on; colorbar; colormap turbo;

% --- Final surface on grid (3D)
figure('Color','w');
surf(xGridMesh, yGridMesh, displ_LSC_grid, 'EdgeColor','none', 'FaceAlpha', 0.85);
hold on;
scatter3(x_sort, y_sort, displ_epc, 8, 'k', 'filled'); % observation points
title('LSC Final Surface (grid) with observations');
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Mean Velocity (mm/yr)');q
colorbar; colormap turbo; grid on; view(3);

% --- Spatial comparison of final fields at observation points
figure('Color','w');
subplot(1,2,1)
scatter(x_sort - mean(x_sort), y_sort - mean(y_sort), 20, displ_coll - displ_LSC_final, 'filled');
title('Diff: (OLS+Coll) - (LSC final)'); grid on; colorbar; colormap turbo;
xlabel('X (m, zero-mean)'); ylabel('Y (m, zero-mean)');

subplot(1,2,2)
scatter(x_sort - mean(x_sort), y_sort - mean(y_sort), 20, poly_spl - displ_LSC_final, 'filled');
title('Diff: (Spline) - (LSC final)'); grid on; colorbar; colormap turbo;
xlabel('X (m, zero-mean)'); ylabel('Y (m, zero-mean)');

% --- Optional: Save key results for reporting and reproducibility.
save_LSC = true;
if save_LSC
    save('Etna_LSC_final_results.mat', ...
        'displ_LSC_final','tot_std_LSC','residuals_LSC_final', ...
        'x_sort','y_sort', ...
        'xGridMesh','yGridMesh','displ_LSC_grid', ...
        'poly_LSC2','v_est_LSC','poly_std_LSC2','v_std_LSC', ...
        'corr_LSC_vs_Coll','rms_LSC_vs_Coll','maxabs_LSC_vs_Coll', ...
        'corr_LSC_vs_Spl','rms_LSC_vs_Spl','maxabs_LSC_vs_Spl');
    fprintf('\nSaved LSC final products to Etna_LSC_final_results.mat\n');
end
