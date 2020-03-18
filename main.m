% main.m
clear;
warning('off','all')
addpath('./libs/');
addpath('./lldistkm/');

tic
%% pre-process
% get the mean, var and count of each type of data
f_list = dir('./data/*Primary*.csv'); % use all primary data
data_list = preprocess(f_list);

%% obtain covariance matrix and correlation matrix of certain data types
[cov_mat_pm1, corr_mat_pm1] = cov_matrix(f_list, 'pm1', max(horzcat(data_list(:).pm1_cnt)));
[cov_mat_pm2_5, corr_mat_pm2_5] = cov_matrix(f_list, 'pm2_5', max(horzcat(data_list(:).pm2_5_cnt)));
[cov_mat_pm10, corr_mat_pm10] = cov_matrix(f_list, 'pm10', max(horzcat(data_list(:).pm10_cnt)));
toc

%% fit the RBF kernel
% iteration through each location pair
nPair = (length(f_list) * (length(f_list) - 1)) / 2; % total number of pairs
x = zeros(1, nPair);
y = zeros(1, nPair);
for l1 = 1:length(f_list)
    for l2 = l1+1:length(f_list)
        latlon1=[data_list(l1).lat data_list(l1).lon];
        latlon2=[data_list(l2).lat data_list(l2).lon];
        [d1km, d2km] = lldistkm(latlon1, latlon2);
        x = d1km * 1000; % use d1 distance, convert to meters
        y = cov_mat_pm2_5(l1, l2);
    end
end
f = fit(x,y,'exp1');
plot(f, x, y);

%% plot
lat = horzcat(data_list(:).lat);
lon = horzcat(data_list(:).lon);
sizedata = horzcat(data_list(:).pm2_5_avg);
plot(lat, lon, sizedata);

function plot(lat, lon, sizedata)
    % plot the locations
    figure('Position', [0 0 1000 800]);
    gb = geobubble(lat, lon, sizedata);
    geobasemap streets-light; % set base map style
end
