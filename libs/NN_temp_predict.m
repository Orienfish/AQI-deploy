function [mean_vd] = NN_temp_predict(Xv, Xd, mean_d)
%% Predict the temperature at a new location based on Nearest Neighbor.
%
% Args:
%   Xv: list of locations to predict, [lat lon]
%   Xd: list of locations we observed, [lat lon]
%   mean_d: a vector of the mean of the observations at Xd
%
% Return:
%   mean_vd: a vector of the predicted mean at Xv given observations at Xa

n_V = size(Xv, 1);
n_D = size(Xd, 1);
% a matrix to store the absolute distance between Xv and Xd
dist = zeros(n_V, n_D);
% fill in the matrix
for l1 = 1:n_V
    for l2 = 1:n_D
        dist(l1, l2) = abs(Xv(l1, 1) - Xd(l2, 1)) + ...
            abs(Xv(l1, 2) - Xd(l2, 2));
    end
end
% a column vector containing the minimum value of each row
[min_dist, min_idx] = min(dist, [], 2); 
% fill in the predicted temperature
mean_vd = zeros(n_V, 1);
for l = 1:n_V
    mean_vd(l) = mean_d(min_idx(l));
end
end

