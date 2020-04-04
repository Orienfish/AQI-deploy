function [out] = pSPIEL(Qparams, params)
%% call the pSPIEL provided by SFO toolbox.
%
% Args:
%   Qparams.Xv: list of reference locations to predict, [lat lon]
%   Qparams.cov_vd: cov matrix at Xv given pre-deployment D
%
%   params.n_V: number of reference locations
%   params.m_A: number of sensors to deploy
%   params.K: the fitted RBF kernel function
%   params.logging: logging flag
%
% Return:
%
addpath('./libs/');
addpath('./mlibs/');
addpath('./lldistkm/');
addpath('./gp/');
addpath('./SFO/sfo');

% the ground set for the submodular functions
V_sigma = 1:params.n_V; 
% Mutual information: F_mi(A) = H(V\A) - H(V\A | A)
F_var = sfo_fn_varred(Qparams.cov_vd,V_sigma);
V = V_sigma;

% create the undirected dist matrix
D = zeros(params.n_V); % init a symmetric matrix for connection graph
nodes = Qparams.Xv;
for p = 1:params.n_V
    for q = p+1:params.n_V
    % check the distance to c
    [d1km, d2km] = lldistkm(nodes(p, :), nodes(q, :));
    % update the undirected graph
    D(p, q) = d1km;
    D(q, p) = d1km;
    end
end
D = D + ones(length(Qparams.Xv)); % cost of links + cost of nodes

Q = 0.6*F_var(V); % Quota setting
fprintf("Quota: %f\n", Q);

AP = sfo_pspiel(F_var,V,Q,D); % call pSPIEL
disp(AP);
utility_pspiel = F_var(AP);
[cost_pspiel edges_pspiel steiner_pspiel]= sfo_pspiel_get_cost(AP,D);
disp(sprintf('pSPIEL: Utility = %f, Cost = %f.',utility_pspiel,cost_pspiel));

out.F = utility_pspiel;
out.C = cost_pspiel;
end

