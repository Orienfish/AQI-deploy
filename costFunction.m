function [cost] = costFunction(Qparams, params)
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
%   cost: the cost of the current deployment A
addpath('./mlibs/');
addpath('./lldistkm/');

stack_idx = zeros(params.m_A, 1); % a list of indexes in comm. range as a stack
stack_size = 0;  % init the size of the stack
conn_dist = Inf(params.m_A, 1); % distance from node to sink
commMST = NaN(params.m_A + 1); % init a matrix for connection graph
                               % first n entries are for Xv, the last one 
                               % is for the sink c. A single-direction MST.
predMST = NaN(params.m_A + 1, 1);  % predecessor nodes of the MST
                        
% find the min spanning tree and fill commMST
% get the valid indexes directly connected to the sink
for p = 1:m_A
    % check the distance to c
    [d1km, d2km] = lldistkm(Xa(p, :), params.c);
    if d1km < params.R
        % add the sensor to stack
        stack_size = stack_size + 1;
        stack_idx(stack_size) = p;
        
        conn_dist(p) = d1km; % update distance of the current sensor to sink
        commMST(params.m_A+1, p) = d1km; % update comm. graph
        predMST(p) = params.m_A + 1; % update the predecessor to sink
    end
end
fprintf('init size of stack: %d\n', stack_size);

ptr = 1; % init the ptr for stack
% search in BFS manner on all frontiers
% complexity: O(m_A^2)
while ptr <= stack_size
    % find all nonexplored sensors connected to this one
    cur_idx = stack_idx(ptr);
    for i = 1:m_A        
        % check the distance to the selected sensor in this round
        [d1km, d2km] = lldistkm(Qparams.Xa(cur_idx, :), ...
            Qparams.Xa(i, :));
        if d1km < params.R % connectable
            cur_dist_i = conn_dist(cur_idx) + d1km; % distance to sink
            if isinf(conn_dist(i))
                % first time adding the sensor to stack
                stack_size = stack_size + 1;
                stack_idx(stack_size) = p;
            end
            if cur_dist_i < conn_dist(i)
                % update distance of the current sensor
                conn_dist(i) = cur_dist_i; % update distance to sink
                commMST(cur_idx, i) = d1km; % update comm. graph
                predMST(i) = cur_idx; % update the predecessor
            end
        end
    end
    ptr = ptr + 1; % pop out the first element in the stack
end

connected = ~isinf(predMST); % logical array of connected sensors
if params.logging
    disp('connected status given A:');
    disp(connected);
    disp('commMST:');
    disp(commMST);
    disp('predMST:');
    disp(predMST);
end

% calculate the conditional entropy of all sensors
F = cond_entropy_d(Qparams.Xv, Qparams.cov_vd, Qparams.Xa, Qparams.cov_ad, ...
        params.K);

% calculate the maintenance cost of connected sensors
M = maintain_cost(Qparams.Xa, Qparams.Ta, commMST, predMST, params.logging);

% calculate the penalty of unconnected sensors
P = params.penalty * (params.m_A - sum(connected));

% final cost
cost = params.weights(1) * F + params.weights(2) * M + params(3) * P;
if params.logging
    disp(['cond entropy: ' num2str(F) ' main cost: ' num2str(M) ...
        ' penalty: ' num2str(P)]);
    disp(['total cost: ' num2str(cost)]);
end
end

