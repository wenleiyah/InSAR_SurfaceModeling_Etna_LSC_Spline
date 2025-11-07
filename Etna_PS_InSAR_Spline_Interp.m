clc
clear
close all
warning off

addpath('./Matlab_Functions')
addpath('./SAR_Data')

% ----------------------------- SPLINES - 2D -----------------------------%

%% 1. SAR data 2D - Etna SE
% --- 1.1 Create working directories --- %
% They are required for geoSplinter software
mkdir('data_input')
mkdir('data_output')
mkdir('job')

% --- 1.2 Import SAR data ---%
data = readmatrix('PS_Etna_SE_new.csv');
easting = data(2:end,2);
northing = data(2:end,3);
proj_from = projcrs(3035);
[lat, lon] = projinv(proj_from, easting, northing);

mean_vel = data(2:end,4); % mean_velocity (unit: mm/year)
mean_vel_std = data(2:end,5); % mean_velocity_standard_deviation （mm/year）
displ2D = data(2:end, 6:end);
day_rel = data(1,6:end);

% check the lat/lon
disp([min(lat), max(lat)])
disp([min(lon), max(lon)])

% Convert latitude and longitude in x and y
[x, y] = deg2utm(lat,lon);

%% initial plotting
% mean_velocity geographic scatter plot
figure
geoscatter(lat, lon, 20, mean_vel, 'filled');
colormap(jet);
geobasemap satellite;
colorbar;
title('Mean Velocity of Persistent Scatterers (Etna)');
legend('Mean Velocity (mm/year)', 'Location', 'best');

% Displacement scatter plot at a certain moment (e.g., the 120th epoch)
figure
geoscatter(lat, lon, 20, displ2D(:,120), 'filled'); 
colormap(jet);
geobasemap satellite;
colorbar;
title('Displacement at Epoch 120 ');
legend('Displacement (mm)', 'Location', 'best');

% 3D plot: mean_velocity
figure
plot3(lat, lon, mean_vel, 'x');
grid on;
xlabel('Latitude');
ylabel('Longitude');
zlabel('Mean Velocity (mm/year)');
title('3D Scatter of Mean Velocity');

% 3D plot:Displacement scatter plot at a certain moment (e.g., the 120th epoch)
figure
plot3(lat, lon, displ2D(:,120), 'x');
grid on;
xlabel('Latitude');
ylabel('Longitude');
zlabel('Displacement (mm)');
title('3D Scatter of Displacement at Epoch 120');
%% Linear regression analysis of individual persistent scatterers (PS) 
% is used to explain the trend of SAR displacement over time, i.e., to 
% evaluate the linear deformation rate (mean velocity) of a given point and 
% its uncertainty.

% Choose the mean velocity as modeling object
displ_epc = mean_vel;

% --- 1.3 Linear regression
displ_PS = displ2D(120,:);   % take a random PS as example

% interpoltion of one persistence scatter
%figure
%plot(day_rel,displ_PS)

% Construct design matrix A, representing the linear model.
Aps = [day_rel' ones(length(day_rel),1)]; 

% Least squares solution, calculate:
% xx_est(1) = slope (velocity estimate, unit mm/day) 
% xx_est(2) = intercept (initial deformation variable) 
% pay attention: the unit of the linear regression is mm/day, the
% mean_velocity we got from egms is with unit mm/year
xx_est = (Aps' * Aps) \ (Aps' * displ_PS');

% Residual vector, reflecting the difference between the predicted value and the actual value at each time point
v_est_ps = displ_PS' - Aps * xx_est;

% Estimate the residual variance (𝜎_0^2), i.e., the mean square error for each observation.
s02_est_ps = v_est_ps' * v_est_ps / (size(Aps,1)-size(Aps,2));

% Calculate the covariance matrix of the slope and intercept: reflects the uncertainty of these estimates.
Cxx_est_ps = s02_est_ps * inv(Aps' * Aps);

% Take the square root of the covariance of the velocity (slope) component → obtain its standard deviation → convert to annual error
% Cxx_est_ps(1,1) is the variance estimate of the slope
% *365*4 indicates that you have a time span of 4 years, so the final result is the total error of the deformation rate estimate for these 4 years
s_est_m = sqrt(Cxx_est_ps(1,1)) * 365; % variance of the whole period, more or less 4 years

%% !Setup the priori sigma
s_teacher = sqrt(mean((mean_vel_std).^2));   % mm/year
% typically we don't average standard deviation, we average the variance, and the compute the square root of variance
s02_apr = (s_teacher * 1.3)^2        % mm²
% s02_apr = 0.3^2;                     % mm²
% relaxing factor, it can be multiplied by 2, larger the variance, larger the variablility of the
% traight lines we admitting, more possible model errors, to enlarge the accepted models
% passing through the data

% The assumed observed error variance is used as the upper limit of the model error for subsequent χ² tests
% when using the mean velocity as the evaluation item, the interpolation
% of the mean_velocity is approximately, the error will be propagated
% through the processing, here we need to add a variance as the Instrument
% Accuracy, to introduce some model error, higher the sigma wider the
% accepted model to pass all the data

%% 2. Polynomial surface modelling 
% From now on we work on the spatial data associated to the chosen epoch.

% --- 2.1 Least Squares
% Design matrix for a 2nd-degree polynomial
A = [ones(size(x)), x, y, x.^2, x.*y, y.^2];


% LS solution
x_est = (A' * A) \ (A' * displ_epc);   % to avoid numerical instability
y_est = A * x_est;
v_est = displ_epc - y_est;
s02_est = v_est' * v_est / (size(A,1) - size(A,2));
Cxx_est = s02_est * inv(A' * A);


% --- 2.2 Significance test on parameters
alpha = 0.05;
t_x_est = x_est ./ sqrt(diag(Cxx_est));
t_lim = tinv(1 - alpha/2, size(A,1) - size(A,2));

fprintf('\n\nt-student test on the significance of the estimated parameters (it. 0): \n')
for i = 1:length(x_est)
    fprintf('Param. %i, t_obs: %.4f (t_lim: %.4f)\n', i, abs(t_x_est(i)), t_lim);
end

%% --- 2.1b  Upgrade to a 3rd-degree polynomial (cubic)
% Add the third-order polynomial terms: x^3, x^2*y, x*y^2, y^3
A3 = [A, x.^3, x.^2.*y, x.*y.^2, y.^3];   % A is the design matrix from previous quadratic model

x_est3 = A3 \ displ_epc;                  % Parameters (10×1)
y_est3 = A3 * x_est3;                     % Fitted values
v_est3 = displ_epc - y_est3;              % Residuals


% Residual variance (unit: (mm/yr)^2)
dof3    = size(A3,1) - size(A3,2);
s02_est3 = (v_est3' * v_est3) / dof3;

%  Covariance of the estimated parameters
Cxx_est3 = s02_est3 * ((A3' * A3) \ eye(size(A3,2)));

% --- 2.2b  t-Student significance test (for each coefficient)
alpha     = 0.05;
SE3       = sqrt(diag(Cxx_est3));
t_x_est3  = x_est3 ./ SE3;
t_lim3    = tinv(1 - alpha/2, dof3);

fprintf('\n[T-test] cubic model (d=3), dof=%d, t_lim=%.4f\n', dof3, t_lim3);
for i = 1:length(x_est3)
    fprintf('beta_%d: t_obs=%.4f\n', i-1, t_x_est3(i));
end

% --- 2.3  Partial F-test versus quadratic surface (test if deg2 → deg3 gives a significant improvement)
% Use the residuals and degrees of freedom from the quadratic model
RSS2 = v_est'  * v_est;                   % Residual sum of squares (quadratic model)
RSS3 = v_est3' * v_est3;                  % Residual sum of squares (cubic model)
p2   = size(A,2);                         % = 6
p3   = size(A3,2);                        % = 10
dfn  = p3 - p2;                           % numerator degrees of freedom
dfd  = dof3;                              % denominator degrees of freedom (use the larger model)
Fval = ((RSS2 - RSS3)/dfn) / (RSS3/dfd);
pval = 1 - fcdf(Fval, dfn, dfd);
fprintf('\n[Partial-F] deg2 -> deg3: F=%.3f, p=%.3g  (df=%d,%d)\n', Fval, pval, dfn, dfd);

% --- 2.4  Optional: Chi-square test on the cubic model (reuse s02_apr)

% s02_apr must be the prior variance ((mm/yr)^2)
chi2_obs3 = s02_est3 / s02_apr;
chi2_low3 = chi2inv(alpha/2,   dof3) / dof3;
chi2_hi3  = chi2inv(1-alpha/2, dof3) / dof3;
fprintf('[Chi^2] cubic: %.4f  (%.4f, %.4f)\n', chi2_obs3, chi2_low3, chi2_hi3);


%%
% --- 2.4 Grid solution
% Solution over a grid
xGrid = min(x):100:max(x);
yGrid = min(y):100:max(y);
[xGridMesh, yGridMesh] = meshgrid(xGrid, yGrid);

% Flatten the grid matrices
xGrid1 = xGridMesh(:);
yGrid1 = yGridMesh(:);

% Design matrix for the grid
A_grid = [ones(size(xGrid1)), xGrid1, yGrid1, xGrid1.^2, xGrid1.*yGrid1, yGrid1.^2];

% Compute the modeled displacements on the grid
y_est_grid = A_grid * x_est;

% Reshape the results back to the grid shape
y_est_grid = reshape(y_est_grid, size(xGridMesh));

% Plot
figure
scatter3(x, y, displ_epc, 10, 'filled')
hold on
surf(xGridMesh, yGridMesh, y_est_grid, 'EdgeColor', 'none', 'FaceAlpha', 0.7)
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Mean Velocity (mm/year)');
title('Polynomial Trend Surface and Observed Data');
colorbar;

% Compute the residuals
res_epc = displ_epc - y_est;


% --- 2.5 TIN
T = delaunay(x, y);
TO = triangulation(T, x, y, y_est);

% Plot
figure
scatter3(x, y, displ_epc, 10, 'filled')
hold on
trisurf(TO, 'EdgeColor', 'none', 'FaceAlpha', 0.7)
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Mean Velocity (mm/year)');
title('TIN Surface (Polynomial)');


% --- 2.6 Test on the model
chi2_obs = s02_est / s02_apr;
chi2_lim2 = chi2inv(1-alpha/2, size(A,1)-size(A,2)) / (size(A,1)-size(A,2));
chi2_lim1 = chi2inv(alpha/2, size(A,1)-size(A,2)) / (size(A,1)-size(A,2));
fprintf('Test on polynomial model, X^2_obs: %.4f (X^2_lim_inf: %.4f, X^2_lim_sup: %.4f)\n', abs(chi2_obs), chi2_lim1, chi2_lim2);
% mean of chi-square is the number of freedom, mean(chi-suqare/n-m) = 1;
% here are actually comparing the chi2_obs with 1?


% -> not satisfied -> splines interpolation

% might be 3 options
% model error: s02_est is not fitting the estimated instrument error
% the esitimated instruments error is wrong
% ...?

% Possible causes:
% (1) Large model error (s02_est not consistent with the data);
% (2) Underestimated noise variance (s02_apr too small);
% (3) Local high-frequency variations not captured by the polynomial order.

%% 3. Understand geoSplinter for 2D splines
% geoSplinter requires for 2D splines that the input file is create with
% data given in an increasing x-coordinate order
% Sort based on coordinates
[~, idx_sort] = sort(x); 
x_sort = x(idx_sort);
y_sort = y(idx_sort);
res_sort = res_epc(idx_sort);

% Plot data with spatial zero-mean
figure
plot3((x-mean(x))/1000, (y-mean(y))/1000, res_epc, 'x')
grid on
xlabel('X (km)');
ylabel('Y (km)');
zlabel('Residual Mean Velocity (mm/year)');
title('Zero-Mean Residuals Distribution');

% Input file creation - residuals
writematrix([x_sort, y_sort, res_sort], ... 
    strcat('.', filesep, 'data_input', filesep, 'Etna.txt'), 'Delimiter', 'space')
% input [x_location, y_location, reasidual mean velocity]
% x,y are the coordinates in UTM projection

% Input file creation - observations
writematrix([x_sort, y_sort, displ_epc], ... 
    strcat('.', filesep, 'data_input', filesep, 'Etna_orig.txt'), 'Delimiter', 'space')

% Find the dimensions' ratio of the area
dim_ratio = (max(y) - min(y)) / (max(x) - min(x));
    % --> The number of splines along y will be dim_ratio the number of 
    % splines along x

% Size of the boundary box
dim_x = max(x) - min(x)
dim_y = max(y) - min(y)

%% First Biliner Splines - experiment
% --- 3.1 Bilinear splines
% Job file creation
data_dim = 2;                        % Dimension of dataset (1D/2D)
type_spl = 1;                        % Type of splines (bilinear/bicubic)
disc_typ = 1;                        % Type of discretization (1delta/2delta)
% order 0 /order 1
file_inp = strcat('Etna.txt');       % Input filename
file_out = strcat('Etna_bil1');      % Output filename
num_obs  = length(res_epc);          % Number of observations
num_col  = 15;                       % Number of y-nodes (= splines)
num_row  = num_col/round(dim_ratio); % Number of x-nodes (= splines)
% defined by observation of the data structure
% figure;plot(x,y,'.') visualize the bounding box of the data
% [x,y] = UTM coordinates (WGS84), and here 1 degree related to 16 meters
% the interpolation is with resolution of 1 meter
% figure;plot3((x-mean(x))/1000,(y-mean(y))/1000,res_epc,'x');grid on;
x_min    = min(x_sort);              % First x-coordinate
x_max    = max(x_sort);              % Last x-coordinate
y_min    = min(y_sort);              % First y-coordinate
y_max    = max(y_sort);              % Last y-coordinate
% lambda   = 0.00022;                   % Regularization parameter (λ) [0 ; 0.01 = (s02_apr/delta^2) / 1km^2)]
% s_teacher = 0.1114 mm/year; s02_apr = (0.1114 * 2.5)^2 ≈ 0.0776 mm²/year²
% figure;plot3((x-mean(x))/1000,(y-mean(y))/1000,res_epc,'x');grid on; 
% to valuate the displacement within 1km, first derivative, delta
% lambda = ( the esitimated sigma / the roughly variation in delta per 1km observed from the dataset) / 1km^2)

    % Firts-tentative lambda
    deltaGrid = mean([dim_x/(num_row-1), dim_y/(num_col-1)]);
    lambda_t = lambdaSplines2D([x_sort, y_sort, res_sort], s02_apr, deltaGrid, type_spl);

lambda   = lambda_t;                 % Regularization parameter (λ)
num_sig  = 8;                        % Number of significant digits

jobFile_analysis2D(data_dim, type_spl, file_inp, file_out, num_obs, ...
num_row, num_col, x_min, x_max, y_min, y_max, lambda, num_sig, disc_typ)

% Job file execution
jobFile_execution(file_out)
[data_bil1, ~, ~, ~] = geoSplinter(strcat('./data_output/', file_out), 'bil');
% - Test χ² on the model
res_bil1 = data_bil1(:,5);

s02_est = sum(res_bil1.^2) / (num_obs-num_col*num_row);
chi2_obs = s02_est / s02_apr;
chi2_lim2 = chi2inv(1-alpha/2, num_obs-num_col*num_row) / (num_obs-num_col*num_row);
chi2_lim1 = chi2inv(alpha/2, num_obs-num_col*num_row) / (num_obs-num_col*num_row);
fprintf('Test on splines model (it. 0), X^2_obs: %.4f (X^2_lim_inf: %.4f, X^2_lim_sup: %.4f)\n', abs(chi2_obs), chi2_lim1, chi2_lim2);

% plot the residual to check
figure
histogram(res_bil1);
xlabel('Residuals (mm/year)');
ylabel('Count');
title('Histogram of Bilinear Spline Residuals (Iteration 0)');

figure
plot3(x_sort, y_sort, res_bil1, 'x');
grid on
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Spline Residual (mm/year)');
title('Bilinear Spline Residuals (Iteration 0)');

%% First Biliner Splines - 30
% --- 3.1 Bilinear splines
% Job file creation
data_dim = 2;                        % Dimension of dataset (1D/2D)
type_spl = 1;                        % Type of splines (bilinear/bicubic)
disc_typ = 1;                        % Type of discretization (1delta/2delta)
% order 0 /order 1
file_inp = strcat('Etna.txt');       % Input filename
file_out = strcat('Etna_bil1');      % Output filename
num_obs  = length(res_epc);          % Number of observations
num_col  = 30;                       % Number of y-nodes (= splines)
num_row  = num_col/round(dim_ratio); % Number of x-nodes (= splines)
% defined by observation of the data structure
% figure;plot(x,y,'.') visualize the bounding box of the data
% [x,y] = UTM coordinates (WGS84), and here 1 degree related to 16 meters
% the interpolation is with resolution of 1 meter
% figure;plot3((x-mean(x))/1000,(y-mean(y))/1000,res_epc,'x');grid on;
x_min    = min(x_sort);              % First x-coordinate
x_max    = max(x_sort);              % Last x-coordinate
y_min    = min(y_sort);              % First y-coordinate
y_max    = max(y_sort);              % Last y-coordinate
% lambda   = 0.00022;                   % Regularization parameter (λ) [0 ; 0.01 = (s02_apr/delta^2) / 1km^2)]
% s_teacher = 0.1114 mm/year; s02_apr = (0.1114 * 2.5)^2 ≈ 0.0776 mm²/year²
% figure;plot3((x-mean(x))/1000,(y-mean(y))/1000,res_epc,'x');grid on; 
% to valuate the displacement within 1km, first derivative, delta
% lambda = ( the esitimated sigma / the roughly variation in delta per 1km observed from the dataset) / 1km^2)

    % Firts-tentative lambda
    deltaGrid = mean([dim_x/(num_row-1), dim_y/(num_col-1)]);
    lambda_t = lambdaSplines2D([x_sort, y_sort, res_sort], s02_apr, deltaGrid, type_spl);

lambda   = lambda_t;                 % Regularization parameter (λ)
num_sig  = 8;                        % Number of significant digits

jobFile_analysis2D(data_dim, type_spl, file_inp, file_out, num_obs, ...
num_row, num_col, x_min, x_max, y_min, y_max, lambda, num_sig, disc_typ)

% Job file execution
jobFile_execution(file_out)
[data_bil1, ~, ~, ~] = geoSplinter(strcat('./data_output/', file_out), 'bil');
% - Test χ² on the model
res_bil1 = data_bil1(:,5);

s02_est = sum(res_bil1.^2) / (num_obs-num_col*num_row);
chi2_obs = s02_est / s02_apr;
chi2_lim2 = chi2inv(1-alpha/2, num_obs-num_col*num_row) / (num_obs-num_col*num_row);
chi2_lim1 = chi2inv(alpha/2, num_obs-num_col*num_row) / (num_obs-num_col*num_row);
fprintf('Test on splines model (it. 0), X^2_obs: %.4f (X^2_lim_inf: %.4f, X^2_lim_sup: %.4f)\n', abs(chi2_obs), chi2_lim1, chi2_lim2);

% plot the residual to check
figure
histogram(res_bil1);
xlabel('Residuals (mm/year)');
ylabel('Count');
title('Histogram of Bilinear Spline Residuals (Iteration 0)');

figure
plot3(x_sort, y_sort, res_bil1, 'x');
grid on
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Spline Residual (mm/year)');
title('Bilinear Spline Residuals (Iteration 1)');

%% TO REMOVE THE OUTLIER

res_epc2 = res_sort;
x_sort2 = x_sort;
y_sort2 = y_sort;

res_epc2(abs(res_bil1)>sqrt(s02_est)*4) = [];
x_sort2(abs(res_bil1)>sqrt(s02_est)*4) = [];
y_sort2(abs(res_bil1)>sqrt(s02_est)*4) = [];

writematrix([x_sort2, y_sort2, res_epc2], ...
    './data_input/Etna_clean1.txt', 'Delimiter', 'space');


%%
% --- 3.1 Bilinear splines - iter.1
% Job file creation
data_dim = 2;                        % Dimension of dataset (1D/2D)
type_spl = 1;                        % Type of splines (bilinear/bicubic)
disc_typ = 1;                        % Type of discretization (1delta/2delta)
% order 0 /order 1
file_inp = strcat('Etna_clean1.txt');       % Input filename
file_out = strcat('Etna_clean1_bil1');      % Output filename
num_obs  = length(res_epc2);          % Number of observations
num_col  = 30;                       % Number of y-nodes (= splines)
num_row = round(num_col / dim_ratio); % Number of x-nodes (= splines)
x_min    = min(x_sort2);              % First x-coordinate
x_max    = max(x_sort2);              % Last x-coordinate
y_min    = min(y_sort2);              % First y-coordinate
y_max    = max(y_sort2);              % Last y-coordinate

    % Firts-tentative lambda
    deltaGrid = mean([dim_x/(num_row-1), dim_y/(num_col-1)]);
    lambda_t = lambdaSplines2D([x_sort2, y_sort2, res_epc2], s02_apr, deltaGrid, type_spl);

lambda   = lambda_t;                 % Regularization parameter (λ)
num_sig  = 8;                        % Number of significant digits

jobFile_analysis2D(data_dim, type_spl, file_inp, file_out, num_obs, ...
num_row, num_col, x_min, x_max, y_min, y_max, lambda, num_sig, disc_typ)

% Job file execution
jobFile_execution(file_out)
[data_bil2, ~, ~, ~] = geoSplinter(strcat('./data_output/', file_out), 'bil');
% - Test χ² on the model
res_bil2 = data_bil2(:,5);

s02_est = sum(res_bil2.^2) / (num_obs-num_col*num_row);
chi2_obs = s02_est / s02_apr;
chi2_lim2 = chi2inv(1-alpha/2, num_obs-num_col*num_row) / (num_obs-num_col*num_row);
chi2_lim1 = chi2inv(alpha/2, num_obs-num_col*num_row) / (num_obs-num_col*num_row);
fprintf('Test on splines model (it. 1), X^2_obs: %.4f (X^2_lim_inf: %.4f, X^2_lim_sup: %.4f)\n', abs(chi2_obs), chi2_lim1, chi2_lim2);

% plot the residual to check
figure
histogram(res_bil2);
xlabel('Residuals (mm/year)');
ylabel('Count');
title('Histogram of Bilinear Spline Residuals (Iteration 1)');

figure
plot3(x_sort2, y_sort2, res_bil2, 'x');
grid on
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Spline Residual (mm/year)');
title('Bilinear Spline Residuals (Iteration 1)');

%% TO REMOVE OUTLIERS – Iteration 2
% Continue outlier removal based on res_bil2 from the first-iter spline fit

res_epc3 = res_epc2;
x_sort3 = x_sort2;
y_sort3 = y_sort2;

res_epc3(abs(res_bil2) > sqrt(s02_est) * 3) = [];
x_sort3(abs(res_bil2) > sqrt(s02_est) * 3) = [];
y_sort3(abs(res_bil2) > sqrt(s02_est) * 3) = [];

writematrix([x_sort3, y_sort3, res_epc3], ...
    './data_input/Etna_clean2.txt', 'Delimiter', 'space');

%% --- 3.1 Bilinear splines - iter.2
% Job file creation
data_dim = 2;                        % Dimension of dataset (1D/2D)
type_spl = 1;                        % Type of splines (bilinear/bicubic)
disc_typ = 1;                        % Type of discretization (1delta/2delta)
% order 0 /order 1
file_inp = strcat('Etna_clean2.txt');       % Input filename
file_out = strcat('Etna_clean2_bil1');      % Output filename
num_obs  = length(res_epc3);          % Number of observations
num_col  = 30;                       % Number of y-nodes (= splines)
num_row = round(num_col / dim_ratio); % Number of x-nodes (= splines)
x_min    = min(x_sort3);              % First x-coordinate
x_max    = max(x_sort3);              % Last x-coordinate
y_min    = min(y_sort3);              % First y-coordinate
y_max    = max(y_sort3);              % Last y-coordinate

%     % Firts-tentative lambda
%     deltaGrid = mean([dim_x/(num_row-1), dim_y/(num_col-1)]                 );
%     lambda_t = lambdaSplines2D([x_sort3, y_sort3, res_epc3], s02_apr, deltaGrid, type_spl);
% 
lambda   = lambda_t;                 % Regularization parameter (λ)
%lambda   = 0.1;                 % Regularization parameter (λ)
num_sig  = 8;                        % Number of significant digits

jobFile_analysis2D(data_dim, type_spl, file_inp, file_out, num_obs, ...
num_row, num_col, x_min, x_max, y_min, y_max, lambda, num_sig, disc_typ)

% Job file execution
jobFile_execution(file_out)
[data_bil3, ~, ~, ~] = geoSplinter(strcat('./data_output/', file_out), 'bil');
% - Test χ² on the model
res_bil3 = data_bil3(:,5);

s02_apr_relaxed=0.32;
s02_est = sum(res_bil3.^2) / (num_obs-num_col*num_row);
chi2_obs = s02_est / s02_apr_relaxed;
chi2_lim2 = chi2inv(1-alpha/2, num_obs-num_col*num_row) / (num_obs-num_col*num_row);
chi2_lim1 = chi2inv(alpha/2, num_obs-num_col*num_row) / (num_obs-num_col*num_row);
fprintf('Test on splines model (it. 2), X^2_obs: %.4f (X^2_lim_inf: %.4f, X^2_lim_sup: %.4f)\n', abs(chi2_obs), chi2_lim1, chi2_lim2);

% plot the residual to check
figure
histogram(res_bil3);
xlabel('Residuals (mm/year)');
ylabel('Count');
title('Histogram of Bilinear Spline Residuals (Iteration 2)');

figure
plot3(x_sort3, y_sort3, res_bil3, 'x');
grid on
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Spline Residual (mm/year)');
title('Bilinear Spline Residuals (Iteration 2)');

% --- Visualize the TIN (Triangulated Irregular Network) of the final residuals ---
% First perform a Delaunay triangulation
T_final = delaunay(x_sort3, y_sort3);
T2_final = triangulation(T_final, x_sort3, y_sort3, data_bil3(:,end-1)); % 用最后一次 bilinear 的结果

figure
scatter3(x_sort3, y_sort3, res_epc3, 10, 'filled')
hold on
trisurf(T2_final, 'EdgeColor', 'none', 'FaceAlpha', 0.7)
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Residual Mean Velocity (mm/year)');
title('TIN Surface after Final Bilinear Spline Residuals');
colorbar;

%% Align the retained points
[~, idx3] = ismember([x_sort3 y_sort3], [x_sort y_sort], 'rows');
best_res_epc = data_bil3(:,end-1);

%% Re-estimate the polynomial trend (after removing outliers)
z_in = displ_epc(idx_sort(idx3));    % Displacement (observations after outlier removal)
x_in = x_sort3;
y_in = y_sort3;

% Quadratic polynomial fit
A = [ones(size(x_in)), x_in, y_in, x_in.^2, x_in.*y_in, y_in.^2];

% --- Numerically stable estimation (avoid normal equations) ---
% Parameters
x_est = A \ z_in;                                 % Estimated polynomial coefficients (6×1)
y_est_fin = A * x_est;                            % Fitted polynomial surface
res_poly = z_in - y_est_fin;                      % Residuals (for subsequent collocation)

% Residual variance (unit: (mm/yr)^2)
n = size(A,1); p = size(A,2);
s02_est = (res_poly' * res_poly) / (n - p);

% Parameter covariance (6×6), avoiding inv()
% Use QR factorization of A
[Q,R] = qr(A,0);
% (A'*A)^{-1} = R \ (R' \ I)
Cxx_est = s02_est * (R \ (R' \ eye(p)));

% Propagate to spatial points: covariance of the trend at each point
Cyy_est = A * Cxx_est * A.';                      % [n × n] covariance matrix

% Final fitted surface (trend + best spline residuals)
displ_fin = y_est_fin + best_res_epc;

%% Visualize the final fitted surface (new trend + spline)
figure
scatter3(x_in, y_in, z_in, 10, 'filled'); hold on
surf(xGridMesh, yGridMesh, ...
     griddata(x_in, y_in, displ_fin, xGridMesh, yGridMesh, 'cubic'), ...
     'EdgeColor','none','FaceAlpha',0.7);
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Mean velocity (mm/yr)');
title('Final Fitted Surface (New Polynomial + Bilinear Spline)');
colorbar;


%% Save the variables for subsequent collocation processing

save('Etna_spline_output_cleaned.mat', ...
    'x_sort3', 'y_sort3', ...
    'z_in', ...               % displacement after outlier removal
    'y_est_fin', ...          % polynomial surface after re-fitting
    'best_res_epc', ...       % spline residuals
    'displ_fin', ...          % final modeled field
    'res_poly', ...           % residuals after new polynomial
    'Cyy_est');               % covariance matrix for LSC
