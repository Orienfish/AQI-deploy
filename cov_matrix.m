function [cov_mat, corr_mat] = cov_matrix(label, max_cnt)
% Calculate the covariance matrix and correlation coefficient matrix.
% Repeated reading is performed. Might be able to improve.
% Question: are those zeros valid data or not?
cov_mat = NaN;
f_list = dir('./data/*Primary*.csv');
mat = NaN(max_cnt, length(f_list));
for f_idx = 1:length(f_list)
    % obtain the list of data files
    f_name = f_list(f_idx).name;
    % fprintf('Processing file %s...\n', f_name);
    % read data
    f_path = append(f_list(f_idx).folder, '/', f_name);
    T = readtable(f_path, 'Delimiter', ',', 'HeaderLines', 0);
        
    switch label
        case 'temp'
            valid_idx = ~isnan(T.Temperature_F);
            mat(1:nnz(valid_idx), f_idx) = T.Temperature_F(valid_idx);
        case 'humid'
            valid_idx = ~isnan(T.Humidity__);
            mat(1:nnz(valid_idx), f_idx) = T.Humidity__(valid_idx);
        case 'pm1'
            valid_idx = ~isnan(T.PM1_0_CF1_ug_m3);
            mat(1:nnz(valid_idx), f_idx) = T.PM1_0_CF1_ug_m3(valid_idx);
        case 'pm2_5'
            valid_idx = ~isnan(T.PM2_5_CF1_ug_m3);
            mat(1:nnz(valid_idx), f_idx) = T.PM2_5_CF1_ug_m3(valid_idx);
        case 'pm10'
            valid_idx = ~isnan(T.PM10_0_CF1_ug_m3);
            mat(1:nnz(valid_idx), f_idx) = T.PM10_0_CF1_ug_m3(valid_idx);
        otherwise
            fprintf('Invalid label!');
            return;
    end
end
cov_mat = cov(mat, 'partialrows');
corr_mat = corrcoef(mat, 'Rows', 'pairwise');
end