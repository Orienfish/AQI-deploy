function [cov_va] = gp_predict(Xv, Xa, K)
%% Predict the mean and variance at locations Xv given Xa.
%  Note that although we plan to place sensors at Xa, we don't know the 
%  observations at Xa yet. Therefore we are not able to predict the mean
%  at Xv from the mean at Xa. 
%  So we only estimate the cov through the fitted kernel function.
%
% Args:
%   Xv: list of locations to predict, [lat lon]
%   Xa: list of locations we are supposed to observe, [lat lon]
%   K: the fitted RBF kernel function
%
% Return:
%   mean_va: a vector of the predicted mean at Xv given Xa
%   cov_va: a matrix of the predicted covariances at Xv given Xa

% calculate Sigma_VA, SigmA_AV, Sigma_VV and Sigma_AA
Sigma_VA = gen_Sigma(Xv, Xa, K);
Sigma_VV = gen_Sigma(Xv, Xv, K) + (1e-10) * (eye(size(Xv, 1)));
Sigma_AA = gen_Sigma(Xa, Xa, K) + (1e-10) * (eye(size(Xa, 1)));

% calculate mean vector and covariance matrix
cov_va = Sigma_VV - (Sigma_VA * Sigma_AA \ Sigma_VA');
end

