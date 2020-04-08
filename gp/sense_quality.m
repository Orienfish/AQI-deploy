function [senseQuality] = sense_quality(Xv, cov_vd, Xa, cov_ad, K)
%% Compute the sensing quality at Xv given Xa, all given pre-deployment Xd.
%  Use the pre-deployment variance data.
%
% Args:
%   Xv: list of reference locations to predict, [lat lon]
%   cov_vd: cov matrix at Xv given pre-deployment D
%   Xa: list of locations we are supposed to observe, [lat lon]
%   cov_ad: cov matrix at Xa given pre-deployment D
%   K: the fitted RBF kernel function
%
% Return:
%   senseQuality: the predicted sensing quality at Xv given Xa

% add noises
cov_vd = cov_vd + (1e-10) * (eye(size(cov_vd, 1)));
cov_ad = cov_ad + (1e-10) * (eye(size(cov_ad, 1)));

% calculate the entropy
H_Xv = 1/2 * log2(det(cov_vd)) + 1/2 * log2((2 * pi * exp(1))^size(cov_vd,1));

% calculate the conditional entropy
% calculate the covariance matrix at Xv given Xa
Sigma_VA = gen_Sigma(Xv, Xa, K);
cov_va = cov_vd - Sigma_VA * (cov_ad \ Sigma_VA');
% n = size(cov_va, 1);
H_XvXa = 1/2 * log2(det(cov_va)) + 1/2 * log2((2 * pi * exp(1))^size(cov_va,1));
%fprintf('H_Xv: %f, H_XvXa: %f\n', H_Xv, H_XvXa);

% calculate the sensing quality
senseQuality = H_Xv - H_XvXa;
end

