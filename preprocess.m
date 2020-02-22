function data_list = preprocess()
% Preprocess the existing data in ./data folder.
% Convert historical readings into average, variance and count.
% A struct with name, location, and statistical parameters are returned.
    % use all primary data
    f_list = dir('./data/*Primary*.csv');
    data_list = repmat( struct( ...
        'name', '', ...
        'lat', 0.0, 'lon', 0.0, ...
        'temp_avg', 0.0, 'temp_var', 0.0, 'temp_cnt', 0, ...
        'humid_avg', 0.0, 'humid_var', 0.0, 'humid_cnt', 0, ...
        'pm1_avg', 0.0, 'pm1_var', 0.0, 'pm1_cnt', 0, ...
        'pm2_5_avg', 0.0, 'pm2_5_var', 0.0, 'pm2_5_cnt', 0, ...
        'pm10_avg', 0.0, 'pm10_var', 0.0, 'pm10_cnt', 0), length(f_list));

    for f_idx = 1:length(f_list)
        % obtain the list of data files
        f_name = f_list(f_idx).name;
        % fprintf('Processing file %s...\n', f_name);
        
        % split the file name by '()'
        f_split = strsplit(f_name, {'(', ')'});
        % obtain sensor name
        data_list(f_idx).name = strtrim(char(f_split(1)));
        % obtain latitude and longtitude
        pos = strsplit(char(f_split(4)));
        data_list(f_idx).lat = str2double(pos(1));
        data_list(f_idx).lon = str2double(pos(2));

        % read data
        f_path = append(f_list(f_idx).folder, '/', f_name);
        T = readtable(f_path, 'Delimiter', ',', 'HeaderLines', 0);

        % temperature
        try
            data_list(f_idx).temp_avg = nanmean(T.Temperature_F);
            data_list(f_idx).temp_var = nanvar(T.Temperature_F);
            data_list(f_idx).temp_cnt = nnz(~isnan(T.Temperature_F));
        catch
            fprintf('File %s does not contain temperature\n', f_name);
        end

        % humidity
        try
            data_list(f_idx).humid_avg = nanmean(T.Humidity__);
            data_list(f_idx).humid_var = nanvar(T.Humidity__);
            data_list(f_idx).humid_cnt = nnz(~isnan(T.Humidity__));
        catch
            fprintf('File %s does not contain humidity\n', f_name);
        end

        % PM1
        try
            data_list(f_idx).pm1_avg = nanmean(T.PM1_0_CF1_ug_m3);
            data_list(f_idx).pm1_var = nanvar(T.PM1_0_CF1_ug_m3);
            data_list(f_idx).pm1_cnt = nnz(~isnan(T.PM1_0_CF1_ug_m3));
        catch
            fprintf('File %s does not contain pm1\n', f_name);
        end

        % PM2.5
        try
            data_list(f_idx).pm2_5_avg = nanmean(T.PM2_5_CF1_ug_m3);
            data_list(f_idx).pm2_5_var = nanvar(T.PM2_5_CF1_ug_m3);
            data_list(f_idx).pm2_5_cnt = nnz(~isnan(T.PM2_5_CF1_ug_m3));
        catch
            fprintf('File %s does not contain pm2.5\n', f_name);
        end

        % PM10
        try
            data_list(f_idx).pm10_avg = nanmean(T.PM10_0_CF1_ug_m3);
            data_list(f_idx).pm10_var = nanvar(T.PM10_0_CF1_ug_m3);
            data_list(f_idx).pm10_cnt = nnz(~isnan(T.PM10_0_CF1_ug_m3));
        catch
            fprintf('File %s does not contain pm10\n', f_name);
        end
    end
end