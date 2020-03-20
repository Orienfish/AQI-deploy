function [cov_mat, corr_mat] = cov_matrix(f_list, label, sdate, edate, ...
    interval, cov_mat_sav, corr_mat_sav)
% Calculate the covariance matrix and correlation coefficient matrix.
% Use the data between start date and end date.
% 
% Args:
%   f_list: list of files to import data from
%   label: indicate the type of data to process
%   sdate: the start date in string format 'dd-mm-yyyy HH-MM-SS UTC'
%   edate: the end date in string format 'dd-mm-yyyy HH-MM-SS UTC'
%   interval: the sampling interval in seconds
%   cov_mat_sav: the name of file to save extracted cov mat
%   corr_mat_sav: the name of file to save extracted corr mat
%
% Return:
%   cov_mat: covariance matrix of the data samples of the given label
%   corr_mat: correlation coefficient matrix of the data samples
% 
% Repeated reading is performed. Might be able to improve.
% Question: are those zeros valid data or not? 

% caculate the total length
formatIn = 'yyyy-mm-dd HH:MM:SS UTC';
st_vec = datevec(datenum(sdate, formatIn)); % start time vector
et_vec = datevec(datenum(edate, formatIn)); % end time vector
t_secs = etime(et_vec, st_vec); % calculate the elapsed time between two vectors
n_obs = t_secs / interval + 1; % count of all samples

% reserve a mat to fill in the data of the same type
% col represent variables and rows represent observations
mat = NaN(n_obs, length(f_list));
fprintf('observation mat size %d x %d\n', n_obs, length(f_list));

for f_idx = 1:length(f_list)
    % obtain the list of data files
    f_name = f_list(f_idx).name;
    % fprintf('Processing file %s...\n', f_name);
    % read data
    f_path = append(f_list(f_idx).folder, '/', f_name);
    % we only need one type of data from each file
    T = readtable(f_path, 'Delimiter', ',', 'HeaderLines', 0);
    
    for i = 1:height(T)
        % convert the date to index
        cdate = T.created_at(i);
        curt_vec = datevec(datenum(cdate, formatIn)); % current time vector
        t_secs = etime(curt_vec, st_vec);
        idx = t_secs / interval + 1;
        %fprintf("%d %d\n", t_secs, idx);
        switch label
            % the function can only deal with one label at one time
            case 'temp'
                % fill in the data
                mat(idx, f_idx) = T.Temperature_F(i);
            case 'humid'
                mat(idx, f_idx) = T.Humidity__(i);
            case 'pm1'
                mat(idx, f_idx) = T.PM1_0_CF1_ug_m3(i);
            case 'pm2_5'
                mat(idx, f_idx) = T.PM2_5_CF1_ug_m3(i);
            case 'pm10'
                mat(idx, f_idx) = T.PM10_0_CF1_ug_m3(i);
            otherwise
                fprintf('Invalid label!');
                return;
        end
    end
end
%% cov(A) function: 
% If A is a matrix whose columns represent random 
% variables and whose rows represent observations, C is the covariance 
% matrix with the corresponding column variances along the diagonal.
% 'partialrows' omit rows containing NaN only on a pairwise basis for 
% each two-column covariance calculation
cov_mat = cov(mat, 'partialrows');
writematrix(cov_mat, cov_mat_sav); % write to file

%% corrcoef(A) returns the matrix of correlation coefficients for A, where 
% the columns of A represent random variables and the rows represent observations.
% Use 'pairwise' to compute each two-column correlation coefficient on a 
% pairwise basis. If one of the two columns contains a NaN, that row is omitted.
corr_mat = corrcoef(mat, 'Rows', 'pairwise');
writematrix(corr_mat, corr_mat_sav); % write to file
end