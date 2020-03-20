%% main.m
clear;
warning('off','all')
addpath('./libs/');
addpath('./lldistkm/');

%% settings
% D - pre-deployment
% V - reference locations
% A - deployment plan
m_A = 10;
sdate = '2019-01-01 00:00:00 UTC'; % start date of the dataset
edate = '2020-02-20 23:50:00 UTC'; % end date of the dataset
interval = 60 * 10; % 10 mins

%% pre-process
% get the mean, var and count of each type of data
f_list = dir('./data/*Primary*.csv'); % use all primary data
dataT_sav = './data/dataT.csv';
tic
if exist(dataT_sav, 'file')
    dataT = readtable(dataT_sav);
else
    % if no file exists yet, do pre-process and save it to file
    dataT = preprocess(f_list, dataT_sav);
end
toc
% bubble plot for pre-deployment D on pm2.5
D_lat = vertcat(dataT.lat(:));
D_lon = vertcat(dataT.lon(:));
D = horzcat(D_lat, D_lon);
sizedata = horzcat(dataT.pm2_5_avg(:));
bubbleplot_wsize(D(:, 1), D(:, 2), sizedata);

%% get the region range, set grid unit and obtain V
latUpper = max(D_lat);
latLower = min(D_lat);
lonUpper = max(D_lon);
lonLower = min(D_lon);
gridUnit = 0.05; % adjustable
n_latV = floor((latUpper - latLower) / gridUnit);
n_lonV = floor((lonUpper - lonLower) / gridUnit);
V_lat = linspace(latLower, latUpper, n_latV);
V_lon = linspace(lonLower, lonUpper, n_lonV);

% fill the grid locations into V
n_V = n_latV * n_lonV; % total number of locations in V
V = zeros(n_V, 2);
for i = 1:n_latV
    for j = 1:n_lonV
        V((i - 1) * n_lonV + j, :) = [V_lat(i) V_lon(j)];
    end
end
fprintf('Generate V with size %d x %d \n', n_latV, n_lonV);
% bubbleplot(V(:, 1), V(:, 2));

%% obtain mean vector and covariance matrix and correlation matrix of certain data types
tic
%mean_pm1 = vertcat(dataT.pm1_avg(:));
mean_pm2_5 = vertcat(dataT.pm2_5_avg(:));
%mean_pm10 = vertcat(dataT.pm10_avg(:));
%[cov_mat_pm1, corr_mat_pm1] = cov_matrix(f_list, 'pm1', sdate, edate, interval);
[cov_mat_pm2_5, corr_mat_pm2_5] = cov_matrix(f_list, 'pm2_5', sdate, edate, interval);
%[cov_mat_pm10, corr_mat_pm10] = cov_matrix(f_list, 'pm10', sdate, edate, interval);
toc
figure();
h = heatmap(cov_mat_pm2_5);


%% fit the RBF kernel
% iteration through each location pair
nPair = (length(f_list) * (length(f_list) - 1)) / 2; % total number of pairs
x = zeros(nPair, 1);
y = zeros(nPair, 1);
idx = 1;
for l1 = 1:length(f_list)
    for l2 = l1+1:length(f_list)
        latlon1 = [dataT.lat(l1) dataT.lon(l1)];
        latlon2 = [dataT.lat(l2) dataT.lon(l2)];
        [d1km, d2km] = lldistkm(latlon1, latlon2);
        x(idx) = d1km * 1000; % use d1 distance, convert to meters
        y(idx) = cov_mat_pm2_5(l1, l2);
        idx = idx + 1;
    end
end
% General model Exp1: f(x) = a*exp(b*x)
K = fit(x, y, 'exp1');
figure();
plot(K, x, y);

%% Get the estimated mean and cov at V
[mean_va, cov_va] = gp_predict_knownA(V, D, mean_pm2_5, cov_mat_pm2_5, K);
figure('Position', [0 0 1000 800]);
gb = geodensityplot(V(:, 1), V(:, 2), mean_va); % plot density plot for mean
figure();
h = heatmap(cov_va); % plot heatmap for the cov matrix

%% Generate a random A
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

%% Evaluate uncertainty or conditional entropy
condEntropy = cond_entropy(V, A, K);
fprintf('conditional entropy %f\n', condEntropy);

%% plot functions
function bubbleplot(lat, lon)
    % plot the locations
    figure('Position', [0 0 1000 800]);
    gb = geobubble(lat, lon);
    %geobasemap streets-light; % set base map style
end

function bubbleplot_wsize(lat, lon, sizedata)
    % plot the locations
    figure('Position', [0 0 1000 800]);
    gb = geobubble(lat, lon, sizedata);
    %geobasemap streets-light; % set base map style
end

function bubbleplot_wcolor(lat, lon, sizedata, colordata)
    % plot the locations
    figure('Position', [0 0 1000 800]);
    gb = geobubble(lat, lon, sizedata, colordata);
    %geobasemap streets-light; % set base map style
end
