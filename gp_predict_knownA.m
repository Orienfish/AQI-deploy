function [mean_va, cov_va] = gp_predict_knownA(Xv, Xa, mean_a, cov_a, K)
%% Predict the mean and variance at locations Xv given obervations at Xa.
%  Note that this function assumes Xa are the pre-deployment locations,
%  so we do know the mean and covariance of Xa observations instead of 
%  estimating them through kernel function K.
%
% Args:
%   Xv: list of locations to predict, [lat lon]
%   Xa: list of locations we observed, [lat lon]
%   mean_a, cov_a: a vector and a matrix of the mean and covariance of the
%                  observations at Xa
%   K: the fitted RBF kernel function
%
% Return:
%   mean_va: a vector of the predicted mean at Xv given observations at Xa
%   cov_va: a matrix of the predicted covariances at Xv given observations
%           at Xa

% calculate Sigma_VA, SigmA_AV and Sigma_VV
Sigma_VA = gen_Sigma(Xv, Xa, K);
Sigma_AV = Sigma_VA';
Sigma_VV = gen_Sigma(Xv, Xv, K);

% calculate mean vector and covariance matrix
mean_va = Sigma_VA * inv(cov_a) * mean_a;
cov_va = Sigma_VV - Sigma_VA * inv(cov_a) * Sigma_AV;
end

