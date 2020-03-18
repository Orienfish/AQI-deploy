function [Sigma] = gen_Sigma(locList1, locList2, K)
% Generate the Sigma matrix given two lists of locations and
% the fitted RBF kernel K.
%  
% Args:
%   locList1, locList2: lists of locations [lat lon]
%   K: the fitted RBF kernel function
%
% Return:
%   Sigma_VA: the generated Sigma_VA matrix with |V| rows and |A| cols

% get the number of rows - size of the first dimension 
n_row = size(locList1, 1);
n_col = size(locList2, 1);
Sigma = zeros(n_row, n_col); % init Sigma to  all zeros

% iteration through each location pair
for l1 = 1:n_row
    for l2 = 1:n_col
        latlon1 = locList1(l1, :); % get the l1th location in V
        latlon2 = locList2(l2, :); % get the l2th location in A
        [d1km, d2km] = lldistkm(latlon1, latlon2);
        x = d1km * 1000; % use d1 distance, convert to meters
        Sigma(l1, l2) = K(x); % use the kernel function
    end
end
end

