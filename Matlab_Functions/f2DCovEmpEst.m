function [sigma, eCov, sigmaGrid, eCovF, Cecf, h] = f2DCovEmpEst(S, x, y, dCov, distance, fig)
% function to compute Empirical Covariance function
%
% SYNTAX: [sigma, eCov, sigmaGrid, eCovF, Cecf, h] = f2DCovEmpEst(S, x, y, dCov, distance, fig)
%
% INPUT:
% - S        --> Signal for which the covariance will be estimated
% - x        --> x [m] / long [rad] vector
% - y        --> y [m] / lat [rad] vector
% - dCov     --> step for averaging the point cloud of covariances
% - distance --> kind of distance: 'spherical' or 'cartesian'
% - fig      --> 0 --> no figure,
%                1 --> only averaged covariance,
%                2 --> averaged + sparse points
%
% OUTPUT:
% - sigma     --> sample distance
% - eCov      --> empirical covariance
% - sigmaGrid --> sample distance of empirical covariance function
% - eCovF     --> empirical covariance function sampled at spherical distance sigmaGrid
% - Cecf      --> cofactor matrix of empirical covariance function
% - h         --> handle of the figures (size depending on fig input)

h = -1;  % initialize the handle to the figures

% compute all the spherical distances
if strcmpi(distance, 'spherical')
    % combination of coordinates of all the points
    [L, ~] = meshgrid(y(:), y(:));
    [P, ~] = meshgrid(x(:), x(:));
    % distance
    sigma = atan2(sqrt((cos(P').*sin(L' - L)).^2+(cos(P).*sin(P')-sin(P).*cos(P').*cos(L' - L)).^2),...
                  sin(P).*sin(P')+cos(P).*cos(P').*cos(L' - L)); % better numerical condition
    sigmaScale = 180/pi; % scale for plotting reasons (units from rad to deg)
    clear L P
elseif strcmpi(distance, 'cartesian')
    % combination of coordinates of all the points
    [X, ~] = meshgrid(y(:), y(:));
    [Y, ~] = meshgrid(x(:), x(:));
    % distance
    sigma = sqrt((X-X').^2 + (Y-Y').^2);
    sigmaScale = 1; % scale for plotting reasons (no units change)
    clear X Y
else
    error('Wrong distance definition');
end

% compute the covariance for all the possible couples of points
[Sall, ~] = meshgrid(S(:), S(:));
eCov = Sall .* Sall';
clear S_all

% define the limits of the classes in term of psi values
sigmaMax   = max(sigma(:));      % maximum distance
sigmaClass = (0:dCov:sigmaMax)'; % border of classes

% intitalize output variables
eCovF = zeros(length(sigmaClass),1);     % empirical covariance
Cecf  = zeros(length(sigmaClass),1);     % cofactor matrix diagonal
sigmaGrid = zeros(length(sigmaClass),1); % grid of class centre
% distance equal to 0
eCovF(1) = mean(eCov(sigma==0));
Cecf(1,1) = (1./numel(eCov(sigma==0))).^2;
% other distances
for i = 2:length(sigmaClass) % starting from 2, since the class 1 is the zero distance
    idx = (sigma > sigmaClass(i-1)) & (sigma <= sigmaClass(i)); % filter points inside the class
    eCovF(i)     = mean(eCov(idx));                             % empirical covariance of the class
    sigmaGrid(i) = mean(sigma(idx));                            % mean distance associated to the class
    Cecf(i)      = 1./sum(idx(:)).^2;                           % cofactor matrix (proportional to number of data)
end

if fig >= 1
    h(1) = figure; plot(sigmaGrid.*sigmaScale,eCovF,'.');
    xlim([0 max(sigmaGrid.*sigmaScale)/2]);
    xlabel('distance');
    title('Empirical covariance function');
end
if fig >= 2
    h(2) = figure;
    plot(sigma(:), eCov(:), '.')
    xlabel('distance');
    title('Covariance point cloud');
end
sigma = sigma(:);
eCov = eCov(:);