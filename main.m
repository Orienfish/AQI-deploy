% main.m
warning('off','all')

data_list = preprocess();
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