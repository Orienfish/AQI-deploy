function [Sigma_VA] = gen_Sigma(V, A, K)
% Generate the Sigma_VA matrix given fitted RBF kernel K.
%  
% Args:
%   V: a list of reference locations [lat lon]
%   A: a list of deployment locations [lat lon]
%   K: the fitted RBF kernel function
%
% Return:
%   Sigma_VA: the generated Sigma_VA matrix with |V| rows and |A| cols

% get the number of rows - size of the first dimension 
n_row = size(V, 1);
n_col = size(A, 1);
Sigma_VA = zeros(n_row, n_col); % init Sigma_VA to all zeros

% iteration through each location pair
for l1 = 1:n_row
    for l2 = 1:n_col
        latlon1 = V(l1, :); % get the l1th location in V
        latlon2 = A(l2, :); % get the l2th location in A
        [d1km, d2km] = lldistkm(latlon1, latlon2);
        x = d1km * 1000; % use d1 distance, convert to meters
        Sigma_VA(l1, l2) = K(x); % use the kernel function
    end
end
end

