function [cov_mat, corr_mat] = cov_matrix(f_list, label, max_cnt)
% Calculate the covariance matrix and correlation coefficient matrix.
% 
% Args:
%   f_list: list of files to import data from
%   label: indicate the type of data to process
%   max_cnt: the maximum valid data samples of the specified label
%
% Return:
%   cov_mat: covariance matrix of the data samples of the given label
%   corr_mat: correlation coefficient matrix of the data samples
% 
% Repeated reading is performed. Might be able to improve.
% Question: are those zeros valid data or not?

% reserve a mat to fill in the data of the same type
mat = NaN(max_cnt, length(f_list));
for f_idx = 1:length(f_list)
    % obtain the list of data files
    f_name = f_list(f_idx).name;
    % fprintf('Processing file %s...\n', f_name);
    % read data
    f_path = append(f_list(f_idx).folder, '/', f_name);
    % we only need one type of data from each file
    T = readtable(f_path, 'Delimiter', ',', 'HeaderLines', 0);
        
    switch label
        % the function can only deal with one label at one time
        case 'temp'
            % get the list of indexes of valid data
            valid_idx = ~isnan(T.Temperature_F);
            % copy the valid data to mat
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
%% cov(A) function: If A is a matrix whose columns represent random 
% variables and whose rows represent observations, C is the covariance 
% matrix with the corresponding column variances along the diagonal.
% 'partialrows' omit rows containing NaN only on a pairwise basis for 
% each two-column covariance calculation
cov_mat = cov(mat, 'partialrows');

%% corrcoef(A) returns the matrix of correlation coefficients for A, where 
% the columns of A represent random variables and the rows represent observations.
% Use 'pairwise' to compute each two-column correlation coefficient on a 
% pairwise basis. If one of the two columns contains a NaN, that row is omitted.
corr_mat = corrcoef(mat, 'Rows', 'pairwise');
end