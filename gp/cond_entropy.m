function [condEntropy] = cond_entropy(Xv, Xa, K)
%% Calculate the conditional entropy at Xv given deployment Xa.
%  Totally depends on the kernel function.
%
% Args:
%   Xv: list of reference locations to predict, [lat lon]
%   Xa: list of locations we are supposed to observe, [lat lon]
%   K: the fitted RBF kernel function
%
% Return:
%   condEntropy: the predicted conditional entropy at Xv given Xa

% calculate the covariance matrix at Xv given Xa
cov_va = gp_predict(Xv, Xa, K);
n = size(cov_va, 1);

% calculate the sensing quality or the conditional entropy
condEntropy = 1/2 * log(det(cov_va)) + 1/2 * n * log(2 * pi * exp(1));
end
