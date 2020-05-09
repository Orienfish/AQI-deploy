function dataT = preprocess(f_list, dataT_sav, thres)
% Preprocess the existing data in ./data folder.
% Convert historical readings into average, variance and count.
% A struct with name, location, and statistical parameters are returned.
%
% Args:
%   f_list: list of file names to import data from
%   dataT_sav: the name of file to save extracted dataT
%   thres: a threshold used to filter out outliers
%
% Return:
%   dataT: statistical pattern of the data in a table

% initialize data table dataT
varNames = {'name', 'lat', 'lon', 'obs_cnt', ...
            'temp_avg', 'temp_var', 'temp_cnt', 'temp_max', 'temp_min', ...
            'humid_avg', 'humid_var', 'humid_cnt', 'humid_max', 'humid_min', ...
            'pm1_avg', 'pm1_var', 'pm1_cnt', 'pm1_max', 'pm1_min', ...
            'pm2_5_avg', 'pm2_5_var', 'pm2_5_cnt', 'pm2_5_max', 'pm2_5_min', ...
            'pm10_avg', 'pm10_var', 'pm10_cnt', 'pm10_max', 'pm10_min'};
varTypes = {'string', 'double', 'double', 'int32', ...
            'double', 'double', 'int32', 'double', 'double', ...
            'double', 'double', 'int32', 'double', 'double', ...
            'double', 'double', 'int32', 'double', 'double', ...
            'double', 'double', 'int32', 'double', 'double', ...
            'double', 'double', 'int32', 'double', 'double'};
dataT = table('Size', [length(f_list) length(varNames)], ...
              'VariableTypes', varTypes, 'VariableNames', varNames);

for f_idx = 1:length(f_list)
    % obtain the list of data files
    f_name = f_list(f_idx).name;
    % fprintf('Processing file %s...\n', f_name);
        
    % split the file name by '()'
    f_split = strsplit(f_name, {'(', ')'});
    % obtain sensor name
    dataT.name(f_idx) = strtrim(char(f_split(1)));
    % obtain latitude and longtitude
    pos = strsplit(char(f_split(4)));
    dataT.lat(f_idx) = str2double(pos(1));
    dataT.lon(f_idx) = str2double(pos(2));

    % read data
    f_path = append(f_list(f_idx).folder, '/', f_name);
    T = readtable(f_path, 'Delimiter', ',', 'HeaderLines', 0);
    
    % getthe number of total observations
    dataT.obs_cnt(f_idx) = height(T);

    % temperature
    try
        array = T.Temperature_F(~isnan(T.Temperature_F)); % filter out nan
        array = array(array >= 0); % filter out outliers
        array = array(array < thres); % filter out outliers
        dataT.temp_avg(f_idx) = mean(array);
        dataT.temp_var(f_idx) = var(array);
        dataT.temp_cnt(f_idx) = length(array);
        dataT.temp_max(f_idx) = max(array);
        dataT.temp_min(f_idx) = min(array);
    catch
        fprintf('File %s does not contain temperature\n', f_name);
    end

    % humidity
    try
        array = T.Humidity__(~isnan(T.Humidity__)); % filter out nan
        array = array(array < thres); % filter out outliers
        dataT.humid_avg(f_idx) = mean(array);
        dataT.humid_var(f_idx) = var(array);
        dataT.humid_cnt(f_idx) = length(array);
        dataT.humid_max(f_idx) = max(array);
        dataT.humid_min(f_idx) = min(array);
    catch
        fprintf('File %s does not contain humidity\n', f_name);
    end

    % PM1
    try
        array = T.PM1_0_CF1_ug_m3(~isnan(T.PM1_0_CF1_ug_m3)); % filter out nan
        array = array(array < thres); % filter out outliers
        dataT.pm1_avg(f_idx) = mean(array);
        dataT.pm1_var(f_idx) = var(array);
        dataT.pm1_cnt(f_idx) = length(array);
        dataT.pm1_max(f_idx) = max(array);
        dataT.pm1_min(f_idx) = min(array);
    catch
        fprintf('File %s does not contain pm1\n', f_name)
    end

    % PM2.5
    try
        array = T.PM2_5_CF1_ug_m3(~isnan(T.PM2_5_CF1_ug_m3)); % filter out nan
        array = array(array < thres); % filter out outliers
        dataT.pm2_5_avg(f_idx) = mean(array);
        dataT.pm2_5_var(f_idx) = var(array);
        dataT.pm2_5_cnt(f_idx) = length(array);
        dataT.pm2_5_max(f_idx) = max(array);
        dataT.pm2_5_min(f_idx) = min(array);
    catch
        fprintf('File %s does not contain pm2.5\n', f_name);
    end

    % PM10
    try
        array = T.PM10_0_CF1_ug_m3(~isnan(T.PM10_0_CF1_ug_m3)); % filter out nan
        array = array(array < thres); % filter out outliers
        dataT.pm10_avg(f_idx) = mean(array);
        dataT.pm10_var(f_idx) = var(array);
        dataT.pm10_cnt(f_idx) = length(array);
        dataT.pm10_max(f_idx) = max(array);
        dataT.pm10_min(f_idx) = min(array);
    catch
        fprintf('File %s does not contain pm10\n', f_name);
    end
end
writetable(dataT, dataT_sav);
end