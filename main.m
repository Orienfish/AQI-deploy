% use all the primary data
f_list = dir('./data/*Primary*.csv');
pos_list = cell(length(f_list), 2);
for f_idx = 1:length(f_list)
    % obtain the list of data files
    f_name = f_list(f_idx).name;

    % obtain latitude and longtitude
    % split the file name by '()', extract lat and lon
    f_split = strsplit(f_name, {'(', ')'});
    pos = strsplit(char(f_split(4)));
    pos_list(f_idx, :) = pos;
end
pos_list = str2double(pos_list);

% plot the locations
lat = pos_list(:, 1);
lon = pos_list(:, 2);
figure('Position', [0 0 1000 8000]);
gb = geobubble(lat, lon);
geobasemap streets-light; % set base map style