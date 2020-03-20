function dataT = preprocess(f_list, dataT_sav)
% Preprocess the existing data in ./data folder.
% Convert historical readings into average, variance and count.
% A struct with name, location, and statistical parameters are returned.
%
% Args:
%   f_list: list of file names to import data from
%   dataT_sav: the name of file to save extracted dataT
%
% Return:
%   dataT: statistical pattern of the data in a struct

% define dataT struct
%dataT = repmat( struct( ...
%    'name', '', ...
%    'lat', 0.0, 'lon', 0.0, 'obs_cnt', 0.0, ...
%    'temp_avg', 0.0, 'temp_var', 0.0, 'temp_cnt', 0, ...
%    'humid_avg', 0.0, 'humid_var', 0.0, 'humid_cnt', 0, ...
%    'pm1_avg', 0.0, 'pm1_var', 0.0, 'pm1_cnt', 0, ...
%    'pm2_5_avg', 0.0, 'pm2_5_var', 0.0, 'pm2_5_cnt', 0, ...
%    'pm10_avg', 0.0, 'pm10_var', 0.0, 'pm10_cnt', 0), 1, length(f_list));
% initialize data table dataT
varNames = {'name', 'lat', 'lon', 'obs_cnt', ...
            'temp_avg', 'temp_var', 'temp_cnt', ...
            'humid_avg', 'humid_var', 'humid_cnt', ...
            'pm1_avg', 'pm1_var', 'pm1_cnt', ...
            'pm2_5_avg', 'pm2_5_var', 'pm2_5_cnt', ...
            'pm10_avg', 'pm10_var', 'pm10_cnt'};
varTypes = {'string', 'double', 'double', 'int32', ...
            'double', 'double', 'int32', ...
            'double', 'double', 'int32', ...
            'double', 'double', 'int32', ...
            'double', 'double', 'int32', ...
            'double', 'double', 'int32'};
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
        dataT.temp_avg(f_idx) = nanmean(T.Temperature_F);
        dataT.temp_var(f_idx) = nanvar(T.Temperature_F);
        dataT.temp_cnt(f_idx) = nnz(~isnan(T.Temperature_F));
    catch
        fprintf('File %s does not contain temperature\n', f_name);
    end

    % humidity
    try
        dataT.humid_avg(f_idx) = nanmean(T.Humidity__);
        dataT.humid_var(f_idx) = nanvar(T.Humidity__);
        dataT.humid_cnt(f_idx) = nnz(~isnan(T.Humidity__));
    catch
        fprintf('File %s does not contain humidity\n', f_name);
    end

    % PM1
    try
        dataT.pm1_avg(f_idx) = nanmean(T.PM1_0_CF1_ug_m3);
        dataT.pm1_var(f_idx) = nanvar(T.PM1_0_CF1_ug_m3);
        dataT.pm1_cnt(f_idx) = nnz(~isnan(T.PM1_0_CF1_ug_m3));
    catch
        fprintf('File %s does not contain pm1\n', f_name)
    end

    % PM2.5
    try
        dataT.pm2_5_avg(f_idx) = nanmean(T.PM2_5_CF1_ug_m3);
        dataT.pm2_5_var(f_idx) = nanvar(T.PM2_5_CF1_ug_m3);
        dataT.pm2_5_cnt(f_idx) = nnz(~isnan(T.PM2_5_CF1_ug_m3));
    catch
        fprintf('File %s does not contain pm2.5\n', f_name);
    end

    % PM10
    try
        dataT.pm10_avg(f_idx) = nanmean(T.PM10_0_CF1_ug_m3);
        dataT.pm10_var(f_idx) = nanvar(T.PM10_0_CF1_ug_m3);
        dataT.pm10_cnt(f_idx) = nnz(~isnan(T.PM10_0_CF1_ug_m3));
    catch
        fprintf('File %s does not contain pm10\n', f_name);
    end
end
writetable(dataT, dataT_sav);
end