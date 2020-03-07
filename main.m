% main.m
clear;
warning('off','all')

%tic
%data_list = preprocess();
%[cov_mat_pm1, corr_mat_pm1] = cov_matrix('pm1', max(horzcat(data_list(:).pm1_cnt)));
%[cov_mat_pm2_5, corr_mat_pm2_5] = cov_matrix('pm2_5', max(horzcat(data_list(:).pm2_5_cnt)));
%[cov_mat_pm10, corr_mat_pm10] = cov_matrix('pm10', max(horzcat(data_list(:).pm10_cnt)));
%toc
%lat = horzcat(data_list(:).lat);
%lon = horzcat(data_list(:).lon);
%sizedata = horzcat(data_list(:).pm2_5_avg);
%plot(lat, lon, sizedata);

%function plot(lat, lon, sizedata)
%    % plot the locations
%    figure('Position', [0 0 1000 800]);
%    gb = geobubble(lat, lon, sizedata);
%    geobasemap streets-light; % set base map style
%end
temp = linspace(-5, 40, 10);
batlife = zeros(10);
for i = 1:10
    batlife(i) = bat_lifetime(750, temp(i), 40, 0.1);
end
figure();
plot(temp, batlife);
ylim([17 18.5])