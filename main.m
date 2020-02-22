% main.m
warning('off','all')
data_list = preprocess();

function plot(lat, lon)
    % plot the locations
    %lat = pos_list(:, 1);
    %lon = pos_list(:, 2);
    figure('Position', [0 0 1000 8000]);
    gb = geobubble(lat, lon);
    geobasemap streets-light; % set base map style
end