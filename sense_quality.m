function [senseQuality] = sense_quality(Xv, cov_vd, Xa, cov_ad, K, c, R)
%% Compute the sensing quality at Xv given Xa, all given pre-deployment Xd.
%  Use the pre-deployment variance data.
%
% Args:
%   Xv: list of reference locations to predict, [lat lon]ß
%   cov_vd: cov matrix at Xv given pre-deployment D
%   Xa: list of locations we are supposed to observe, [lat lon]
%   cov_ad: cov matrix at Xa given pre-deployment D
%   K: the fitted RBF kernel function
%   c: position of the sink in [lat los]
%   R: communication range of the sensors in km
%
% Return:
%   senseQuality: the predicted sensing quality at Xv given Xa

% calculate the entropy
H_Xv = 0.5 * log(det(cov_vd)); % the same part: 0.5 * n * log(2 * pi * exp(1));

% calculate the conditional entropy
% calculate the covariance matrix at Xv given Xa
Sigma_VA = gen_Sigma(Xv, Xa, K);
Sigma_AV = Sigma_VA';
cov_ad_inv = inv(cov_ad);
cov_va = cov_vd - Sigma_VA * cov_ad_inv * Sigma_AV;
% n = size(cov_va, 1);
H_XvXa = 0.5 * log(det(cov_va)); % the same part: 0.5 * n * log(2 * pi * exp(1));

% calculate the sensing quality
senseQuality = H_Xv - H_XvXa;
end

