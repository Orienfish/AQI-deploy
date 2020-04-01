function [out] = costFunction(Qparams, params)
%% Cost function to minimize in PSO.
%
% Args:
%   Qparams.Xv: list of reference locations to predict, [lat lon]
%   Qparams.cov_vd: cov matrix at Xv given pre-deployment D
%   Qparams.Xa: list of locations we are supposed to observe, [lat lon]
%   Qparams.Ta: average temperature estimation at Xa in Celsius
%   Qparams.cov_ad: cov matrix at Xa given pre-deployment D
%   params.m_A: number of sensors to deploy
%   params.Cm: maintenance cost budget
%   params.K: the fitted RBF kernel function
%   params.c: position of the sink in [lat lon]
%   params.R: communication range of the sensors in km
%   params.weights: 1x3 vectors for the weights
%   params.penalty: penalty for non-connected nodes
%   params.logging: logging flag
%
% Return:
%   out.cost: the cost of the current deployment A
%   out.F: sensing quality of the solution
%   out.M: maintenance cost of the solution
%   out.G: the feasible connection graph of the current deployment A
%   out.pred: the predecessor in MST of A
addpath('./mlibs/');
addpath('./lldistkm/');
addpath('./gp/');

G = zeros(params.m_A + 1); % init a symmetric matrix for connection graph

% create the undirected graph and fill the matrix
nodes = vertcat(Qparams.Xa, params.c);
for p = 1:params.m_A+1
    for q = p+1:params.m_A+1
    % check the distance to c
    [d1km, d2km] = lldistkm(nodes(p, :), nodes(q, :));
    if d1km < params.R
        % update the undirected graph
        G(p, q) = d1km;
        G(q, p) = d1km;
    end
    end
end

% find the minimal spanning tree in graph, with the sink as the root
[Tree, pred] = graphminspantree(sparse(G), params.m_A+1);
connected = ~isnan(pred); % a logical array of connected sensors

if params.logging
    disp('G:');
    disp(G);
    disp('Tree:');
    disp(Tree);
    disp('pred:');
    disp(pred);
    disp('connected:');
    disp(connected);
end

% calculate the sensing quality of all sensors
F = sense_quality(Qparams.Xv, Qparams.cov_vd, Qparams.Xa, Qparams.cov_ad, ...
    params.K);
F = real(F); % take the real part

% calculate the maintenance cost of connected sensors
M = maintain_cost(Qparams.Xa, Qparams.Ta, connected, G, pred, ...
    params.logging);

% calculate the penalty of unconnected sensors
P = params.penalty * (sum(~connected));

% final cost
cost = params.weights(1) * (-F) + params.weights(2) * max(M.C - params.Cm, 0) ... 
    + params.weights(3) * P;
if params.logging
    fprintf('sensing quality: %f main cost: %f penalty: %f\n', F, M.C, P);
    fprintf('total cost: %f\n', cost);
end

% fill in the output results
out.cost = cost;
out.F = F;
out.M = M.C;
out.G = G;
out.pred = pred;
end

