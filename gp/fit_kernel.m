function [K] = fit_kernel(lat, lon, cov_mat, type)
% Fit the RBF kernel.
%
% Args:
%   lat: list of latitude of locations, size (n, 1)
%   lon: list of longitude of locations, size (n, 1)
%   cov_mat: the covariance matrix of the values at these locations, size (n, n)
%   type: data type string
%
% Return:
%   K: the fitted kernel function

% iteration through each location pair
n = length(lat); % total number of locations
nPair = (n * (n - 1)) / 2; % total number of pairs
x = zeros(nPair, 1);
y = zeros(nPair, 1);
idx = 1;
for l1 = 1:n
    for l2 = l1+1:n
        % omit nan and negative covariances
        if ~isnan(cov_mat(l1, l2)) && cov_mat(l1, l2) > 0
            latlon1 = [lat(l1) lon(l1)];
            latlon2 = [lat(l2) lon(l2)];
            [d1km, d2km] = lldistkm(latlon1, latlon2);
            x(idx) = d1km * 1000; % use d1 distance, convert to meters
            y(idx) = cov_mat(l1, l2);
            idx = idx + 1;
        end
    end
end
% General model Exp1: f(x) = a*exp(b*x)
K = fit(x, y, 'exp1');
disp(type); K
figure();
plot(K, x, y);
title(type);
end

