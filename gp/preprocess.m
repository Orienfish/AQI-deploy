function dataT = preprocess(f_list, dataT_sav, thres, N_bin)
% Preprocess the existing data in ./data folder.
% Convert historical readings into average, variance and count.
% A struct with name, location, and statistical parameters are returned.
%
% Args:
%   f_list: list of file names to import data from
%   dataT_sav: the name of file to save extracted dataT
%   thres: a threshold used to filter out outliers
%   N_bin: number of bins to approximate ambient temperature distribution
%
% Return:
%   dataT: statistical pattern of the data in a table
addpath('../libs/');
addpath('../mlibs/');

% initialize data table dataT
varNames = {'name', 'lat', 'lon', 'obs_cnt', ...
            'temp_avg', 'temp_var', 'temp_cnt', 'temp_max', 'temp_min', ...
            'humid_avg', 'humid_var', 'humid_cnt', 'humid_max', 'humid_min', ...
            'pm1_avg', 'pm1_var', 'pm1_cnt', 'pm1_max', 'pm1_min', ...
            'pm2_5_avg', 'pm2_5_var', 'pm2_5_cnt', 'pm2_5_max', 'pm2_5_min', ...
            'pm10_avg', 'pm10_var', 'pm10_cnt', 'pm10_max', 'pm10_min'};
varTypes = {'string', 'double', 'double', 'int32', ...
            'double', 'double', 'double', 'double', 'double', ...
            'double', 'double', 'int32', 'double', 'double', ...
            'double', 'double', 'int32', 'double', 'double', ...
            'double', 'double', 'int32', 'double', 'double', ...
            'double', 'double', 'int32', 'double', 'double'};
dataT = table('Size', [length(f_list) length(varNames)], ...
              'VariableTypes', varTypes, 'VariableNames', varNames);
figure('Position', [0 0 500 375]);

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
    if sum(strcmp('Temperature_F', T.Properties.VariableNames)) > 0
        array = T.Temperature_F(~isnan(T.Temperature_F)); % filter out nan
        array = array(array >= 0); % filter out outliers
        array = array(array < thres); % filter out outliers
        dataT.temp_avg(f_idx) = mean(array);
        dataT.temp_var(f_idx) = var(array);
        dataT.temp_cnt(f_idx) = length(array);
        dataT.temp_max(f_idx) = max(array);
        dataT.temp_min(f_idx) = min(array);
        
        array = fah2cel(array);
        if mod(f_idx, 12) == 0
            [curCounts, curEdges] = histcounts(array, N_bin);
            curCenters = 0.5 * (curEdges(1:N_bin) + curEdges(2:N_bin+1));
            curCounts = curCounts ./ dataT.temp_cnt(f_idx);
            h = cdfplot(array);
            set(h, 'LineWidth', 1.5);
            hold on;
            
            % calculate the maintenance cost at this position
            params.Pto = 0.52;       % 520mW transmission power baseline
            params.Btx = 125;        % 1kbps = 125B/s bandwidth
            params.Brx = 125;
            params.dtx = 10;         % 10km maximal comm. range
            params.Ltx = 100;        % 100B packet length
            params.Lrx = 100;           
            params.Prx = 0.2;        % 200mW receiving power
            params.Psen = 0.2;       % 200mW senisng power
            params.tsen = 0.3;       % 300ms sensing time
            params.T = 10;           % 10s sampling frequency
            params.f = 300e6;        % 300MHz clock frequency
            params.Vdd = 3.3;        % 3.3v supply voltage
            params.Iref = 50;        % 50mA reference current draw
            params.Pref = params.Vdd * params.Iref / 1000; % reference power in W
            params.nbins = 10;       % number of bins to deal with temperature variation
            
            % settings for battery
            params.cap_bat = 2000;   % initial battery capacity in mAh
            params.dt_bat_h = 1;     % time resolution of battery in hours
            params.c_bat = 10;       % cost to replace battery

            % setting for circuit
            params.c_node = 100;     % cost to replace node
            
            [nodeC, batlife_ratio, cirlife_ratio] = maintain_node(params, ...
                curCenters, curCounts);
            fprintf('maintain cost: %f batlife: %f cirlife: %f\n', ...
                nodeC, batlife_ratio, cirlife_ratio);
        end
    else
        fprintf('File %s does not contain temperature\n', f_name);
    end

    % humidity
    if sum(strcmp('Humidity__', T.Properties.VariableNames)) > 0
        array = T.Humidity__(~isnan(T.Humidity__)); % filter out nan
        array = array(array < thres); % filter out outliers
        dataT.humid_avg(f_idx) = mean(array);
        dataT.humid_var(f_idx) = var(array);
        dataT.humid_cnt(f_idx) = length(array);
        dataT.humid_max(f_idx) = max(array);
        dataT.humid_min(f_idx) = min(array);
    else
        fprintf('File %s does not contain humidity\n', f_name);
    end

    % PM1
    if sum(strcmp('PM1_0_CF1_ug_m3', T.Properties.VariableNames)) > 0
        array = T.PM1_0_CF1_ug_m3(~isnan(T.PM1_0_CF1_ug_m3)); % filter out nan
        array = array(array < thres); % filter out outliers
        dataT.pm1_avg(f_idx) = mean(array);
        dataT.pm1_var(f_idx) = var(array);
        dataT.pm1_cnt(f_idx) = length(array);
        dataT.pm1_max(f_idx) = max(array);
        dataT.pm1_min(f_idx) = min(array);
    else
        fprintf('File %s does not contain pm1\n', f_name)
    end

    % PM2.5
    if sum(strcmp('PM2_5_CF1_ug_m3', T.Properties.VariableNames)) > 0
        array = T.PM2_5_CF1_ug_m3(~isnan(T.PM2_5_CF1_ug_m3)); % filter out nan
        array = array(array < thres); % filter out outliers
        dataT.pm2_5_avg(f_idx) = mean(array);
        dataT.pm2_5_var(f_idx) = var(array);
        dataT.pm2_5_cnt(f_idx) = length(array);
        dataT.pm2_5_max(f_idx) = max(array);
        dataT.pm2_5_min(f_idx) = min(array);
    else
        fprintf('File %s does not contain pm2.5\n', f_name);
    end

    % PM10
    if sum(strcmp('PM10_0_CF1_ug_m3', T.Properties.VariableNames)) > 0
        array = T.PM10_0_CF1_ug_m3(~isnan(T.PM10_0_CF1_ug_m3)); % filter out nan
        array = array(array < thres); % filter out outliers
        dataT.pm10_avg(f_idx) = mean(array);
        dataT.pm10_var(f_idx) = var(array);
        dataT.pm10_cnt(f_idx) = length(array);
        dataT.pm10_max(f_idx) = max(array);
        dataT.pm10_min(f_idx) = min(array);
    else
        fprintf('File %s does not contain pm10\n', f_name);
    end
end
ax = gca; ax.FontSize=16;
xlim([0, 40]);
xlabel('Temperature (°C)'); ylabel('Cumulative Probability');
writetable(dataT, dataT_sav);
end