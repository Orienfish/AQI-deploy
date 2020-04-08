function [out] = pSPIEL(Qparams, params)
%% call the pSPIEL provided by SFO toolbox.
%
% Args:
%   Qparams.Xv: list of reference locations to predict, [lat lon]
%   Qparams.cov_vd: cov matrix at Xv given pre-deployment D
%   Qparams.mean_temp_d: mean temperature at D
%   Qparams.Xd: list of predeployment locations
%   Qparams.mean_d: mean value at D
%   Qparams.cov_d: cov matrix at D
%   Qparams.Xa: list of locations we are supposed to observe, [lat lon]
%   Qparams.Ta: average temperature estimation at Xa in Celsius
%   Qparams.cov_ad: cov matrix at Xa given pre-deployment D
%   Qparams.mean_temp_d: mean temperature at D
%   Qparams.cov_temp_d: cov matrix of temperature at D
%
%   params.n_V: number of reference locations
%   params.m_A: number of sensors to deploy
%   params.Q: sensing quality quota
%   params.K: the fitted RBF kernel function
%   parmas.K_temp: the fitted RBF kernel function for temperature
%   params.c: position of the sink in [lat lon]
%   params.R: communication range of the sensors in km
%   params.penalty: penalty for non-connected nodes
%   params.logging: logging flag
%
% Return:
%   out.F: sensing quality
%   out.M: maintenance cost
%   out.Position: select locations
%   out.pred: the predecessors of each node
%   out.connected: 0 or 1 indicating which node is selected
addpath('./libs/');
addpath('./mlibs/');
addpath('./lldistkm/');
addpath('./gp/');
addpath('./SFO/sfo');

% the ground set for the submodular functions
V_sigma = 1:params.n_V;
V = V_sigma;

% Mutual information: F_mi(A) = H(V\A) - H(V\A | A)
F_mi = sfo_fn_mi(Qparams.cov_vd,V_sigma);
% Conditional entropy H(V\A | A)
F_ig = sfo_fn_infogain(Qparams.cov_vd,V_sigma,.01);
% Variance reduction: F_var(A) = Var(V)-Var(V | A)
F_var = sfo_fn_varred(Qparams.cov_vd,V_sigma);

% choose one of the criterion and set Quota
% need manual tuning here
choose = 'mi';
if strcmp(choose, 'mi')
    F = F_mi; 
    Q = params.Q;
elseif strcmp(choose, 'ig')
    F = F_ig;
    Q = 0.8 * F(V);
elseif strcmp(choose, 'var')
    F = F_var;
    Q = 0.6 * F(V);
else
    error('wrong criterion selection!');
end
%fprintf("Choose %s, Quota: %f\n", choose, Q);

% create the undirected dist matrix
D = zeros(params.n_V); % init a symmetric matrix for connection graph
nodes = Qparams.Xv;
for p = 1:params.n_V
    for q = p+1:params.n_V
    % check the distance to c
    [d1km, d2km] = lldistkm(nodes(p, :), nodes(q, :));
    if d1km > params.R
        D(p, q) = d1km + params.penalty;
        D(q, p) = d1km + params.penalty;
    else
        % update the undirected graph
        D(p, q) = d1km;
        D(q, p) = d1km;
    end
    end
end
D = D + ones(length(Qparams.Xv)); % cost of links + cost of nodes

[AP, E, res] = sfo_pspiel(F,V,Q,D); % call pSPIEL
%disp(AP);
utility_pspiel = F(AP);
[cost_pspiel edges_pspiel steiner_pspiel]= sfo_pspiel_get_cost(AP, D);
disp(sprintf('pSPIEL: Utility = %f, Cost = %f.',utility_pspiel,cost_pspiel));

% generate the node selection by pSPIEL
select = zeros(params.n_V, 1);
select(AP) = 1;
select = logical(select);

% reconstruct the connection graph
G = zeros(params.n_V + 1); % init a symmetric matrix for connection graph

% fill in D from E
for i=1:size(E, 1)
    p = E(i, 1);
    q = E(i, 2);
    [d1km, d2km] = lldistkm(nodes(p, :), nodes(q, :));
    if d1km < params.R
        G(p, q) = d1km;
        G(q, p) = d1km;
    end
end

% find all nodes adjacent to the sink
for i=1:params.n_V
    if select(i) == 1 % selected node
        [d1km, d2km] = lldistkm(nodes(i, :), params.c);
        if d1km < params.R
            G(i, params.n_V+1) = d1km;
            G(params.n_V+1, i) = d1km;
        end
    end
end
[Tree, pred] = graphminspantree(sparse(G), params.n_V+1);
connected = logical(~isnan(pred(1:params.n_V))); % a logical array of connected sensors
Xa = Qparams.Xv(connected, :);
Xa_remain = Qparams.Xv(~connected, :);
cov_Xa = Qparams.cov_vd(connected, connected);
cov_Xa_remain = Qparams.cov_vd(~connected, ~connected);

% predict the ambient temperature at Xv in Celsius
[temp_mean_vd, temp_cov_vd] = gp_predict_knownD( ...
    Qparams.Xv, Qparams.Xd, Qparams.mean_temp_d, Qparams.cov_temp_d, ...
    params.K_temp);
temp_mean_vd = temp_mean_vd / 4 + 180; % weird fix
Tv = fah2cel(temp_mean_vd);  % convert to Celsius

% update the sensing quality and maintenance cost
out.F = sense_quality(Xa_remain, cov_Xa_remain, Xa, cov_Xa, params.K);
out.M = maintain_cost(Qparams.Xv, Tv, connected, G, pred, false); % one-to-one match
out.Position = Qparams.Xv;
out.pred = pred;
out.connected = connected;
end

