% main.m
clear;
warning('off','all')
addpath('./libs/');
addpath('./lldistkm/');

%% pre-process
% get the mean, var and count of each type of data
f_list = dir('./data/*Primary*.csv'); % use all primary data
tic
data_list = preprocess(f_list); % pre-process
toc
% bubble plot
lat = horzcat(data_list(:).lat);
lon = horzcat(data_list(:).lon);
sizedata = horzcat(data_list(:).pm2_5_avg);
bubbleplot_wsize(lat, lon, sizedata);

%% get the region range, set grid unit and obtain V
latUpper = max(lat);
latLower = min(lat);
lonUpper = max(lon);
lonLower = min(lon);
gridUnit = 0.05; % adjustable
n_latV = floor((latUpper - latLower) / gridUnit);
n_lonV = floor((lonUpper - lonLower) / gridUnit);
latV = linspace(latLower, latUpper, n_latV);
lonV = linspace(lonLower, lonUpper, n_lonV);

% fill the grid locations into V
n_V = n_latV * n_lonV; % total number of locations in V
V = zeros(n_V, 2);
for i = 1:n_latV
    for j = 1:n_lonV
        V((i - 1) * n_lonV + j, :) = [latV(i) lonV(j)];
    end
end
fprintf('Generate V with size %d x %d \n', n_latV, n_lonV);
% bubbleplot(V(:, 1), V(:, 2));

%% obtain covariance matrix and correlation matrix of certain data types
tic
[cov_mat_pm1, corr_mat_pm1] = cov_matrix(f_list, 'pm1', max(horzcat(data_list(:).pm1_cnt)));
[cov_mat_pm2_5, corr_mat_pm2_5] = cov_matrix(f_list, 'pm2_5', max(horzcat(data_list(:).pm2_5_cnt)));
[cov_mat_pm10, corr_mat_pm10] = cov_matrix(f_list, 'pm10', max(horzcat(data_list(:).pm10_cnt)));
toc

%% fit the RBF kernel
% iteration through each location pair
nPair = (length(f_list) * (length(f_list) - 1)) / 2; % total number of pairs
x = zeros(nPair, 1);
y = zeros(nPair, 1);
idx = 1;
for l1 = 1:length(f_list)
    for l2 = l1+1:length(f_list)
        latlon1 = [data_list(l1).lat data_list(l1).lon];
        latlon2 = [data_list(l2).lat data_list(l2).lon];
        [d1km, d2km] = lldistkm(latlon1, latlon2);
        x(idx) = d1km * 1000; % use d1 distance, convert to meters
        y(idx) = cov_mat_pm2_5(l1, l2);
        idx = idx + 1;
    end
end
% General model Exp1: f(x) = a*exp(b*x)
K = fit(x, y, 'exp1');
% plot(K, x, y);

%% Generate a random A
m_A = 10;
A = zeros(m_A, 2);
pd_lat = makedist('Uniform', 'lower', latLower, 'upper', latUpper);
pd_lon = makedist('Uniform', 'lower', lonLower, 'upper', lonUpper);
for i = 1:m_A
    A(i, :) = [random(pd_lat) random(pd_lon)];
end
% plot V and A on one map
sizedata = ones(m_A + n_V, 1); % uniform size
colordata = categorical(vertcat(zeros(m_A, 1), ones(n_V, 1))); % two categories
bubbleplot_wcolor(vertcat(A(:, 1), V(:, 1)), vertcat(A(:, 2), V(:, 2)), ...
    sizedata, colordata);

%% plot functions
function bubbleplot(lat, lon)
    % plot the locations
    figure('Position', [0 0 1000 800]);
    gb = geobubble(lat, lon);
    geobasemap streets-light; % set base map style
end

function bubbleplot_wsize(lat, lon, sizedata)
    % plot the locations
    figure('Position', [0 0 1000 800]);
    gb = geobubble(lat, lon, sizedata);
    geobasemap streets-light; % set base map style
end

function bubbleplot_wcolor(lat, lon, sizedata, colordata)
    % plot the locations
    figure('Position', [0 0 1000 800]);
    gb = geobubble(lat, lon, sizedata, colordata);
    geobasemap streets-light; % set base map style
end
