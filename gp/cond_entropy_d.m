function [condEntropy] = cond_entropy_d(Xv, cov_vd, Xa, cov_ad, K)
%% Calculate the conditional entropy at Xv given deployment Xa.
%  Use the pre-deployment variance data.
%
% Args:
%   Xv: list of reference locations to predict, [lat lon]
%   cov_vd: cov matrix at Xv given pre-deployment Dß
%   Xa: list of locations we are supposed to observe, [lat lon]
%   cov_ad: cov matrix at Xa given pre-deployment D
%   K: the fitted RBF kernel function
%
% Return:
%   condEntropy: the predicted conditional entropy at Xv given Xa

% calculate the covariance matrix at Xv given Xa
Sigma_VA = gen_Sigma(Xv, Xa, K);
Sigma_AV = Sigma_VA';
cov_ad_inv = inv(cov_ad);
cov_va = cov_vd - Sigma_VA * cov_ad_inv * Sigma_AV;

n = size(cov_va, 1);

% calculate the sensing quality or the conditional entropy
condEntropy = 0.5 * log(det(cov_va)) + 0.5 * n * log(2 * pi * exp(1));
end

