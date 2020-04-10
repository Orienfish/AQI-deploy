function [mat] = fix_matrix(locList, mat, K)
%% Fix the covariance matrix with NaN and zero values
%
% Args:
%   locList: the location list
%   mat: the covariance matrix to be fixed
%   K: the generated kernel for the current dataset
for l1 = 1:size(mat, 1)
    for l2 = l1+1:size(mat, 2)
        if isnan(mat(l1, l2)) || mat(l1, l2) < 0
            latlon1 = locList(l1, :); % get the l1th location in V
            latlon2 = locList(l2, :); % get the l2th location in A
            [d1km, d2km] = lldistkm(latlon1, latlon2);
            x = d1km * 1000; % use d1 distance, convert to meters
            mat(l1, l2) = K(x); % use the kernel function
            mat(l2, l1) = K(x);
        end
    end
end
end

