%% To Select the Time Period before the earthquake

% Read relative days (time steps)
% day_rel = data(1, 6:end);   % [0, 6, 12, ..., 1806]
% Start time (knowing that 0 corresponds to 2019-01-13)
% t0 = datetime('2019-01-13');
% Constructing the real time vector
% date_vec = t0 + days(day_rel);  % [2019-01-13, 2019-01-19, ..., ...]
% Find the data index before the earthquake
% idx_pre_eq = date_vec < datetime('2021-09-21');
% Retrieve displacement data before an earthquake
% displ2D_pre = displ2D(:, idx_pre_eq);
% day_rel_pre = day_rel(idx_pre_eq);  % The “relative number of days” before an earthquake
