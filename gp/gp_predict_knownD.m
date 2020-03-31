function [mean_vd, cov_vd] = gp_predict_knownD(Xv, Xd, mean_d, cov_d, K)
%% Predict the mean and variance at locations Xv given obervations at Xd.
%  Note that this function assumes Xd are the pre-deployment locations,
%  so we do know the mean and covariance of observations at Xd.
%
% Args:
%   Xv: list of locations to predict, [lat lon]
%   Xd: list of locations we observed, [lat lon]
%   mean_d, cov_d: a vector and a matrix of the mean and covariance of the
%                  observations at Xd
%   K: the fitted RBF kernel function
%
% Return:
%   mean_va: a vector of the predicted mean at Xv given observations at Xa
%   cov_va: a matrix of the predicted covariances at Xv given observations
%           at Xa

% calculate Sigma_VD, Sigma_DV and Sigma_VV
Sigma_VD = gen_Sigma(Xv, Xd, K);
Sigma_DV = Sigma_VD';
Sigma_VV = gen_Sigma(Xv, Xv, K);

% calculate mean vector and covariance matrix
mean_vd = Sigma_VD * inv(cov_d) * mean_d;
cov_vd = Sigma_VV - Sigma_VD * inv(cov_d) * Sigma_DV;
end

