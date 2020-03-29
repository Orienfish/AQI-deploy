function [Xa, commMST] = IDSQ(Xv, cov_vd, Tv, params)
%% Greedy heuristic IDSQ to place sensors on a subset of locations.
%
% Args:
%   Xv: a list of candidate locations to choose from
%   cov_vd: cov matrix at Xv given pre-deployment D
%   Tv: average temperature estimation at Xv in Celsius
%   params.m_A: number of sensors to deploy
%   params.Cm: maintenance cost budget
%   params.K: the fitted RBF kernel function
%   params.c: position of the sink in [lat lon]
%   params.R: communication range of the sensors in km
%   params.alpha: the weight factor in IDSQ
%
% Return:
%   Xa: a list of solution locations to place sensors
%   commMST: the generated communication graph
addpath('./libs/');
addpath('./mlibs/');
addpath('./lldistkm/');
n_V = size(Xv, 1); % number of reference locations Xv
valid_idx = zeros(n_V, 1); % a list of valid indexes in comm. range
Xa_idx = zeros(n_V, 1); % init idx list for Xa to all zero
lastF = 0.0; % sensing quality in last round of greedy selection
lastM = 0.0; % maintenance cost in last round of greedy selection
commMST = NaN(n_V + 1); % init a matrix for connection graph
                        % first n entries are for Xv, the last one is for 
                        % the sink c. A single-direction MST.
logging = false;

% get the valid indexes directly connected to the sink
for p = 1:length(valid_idx)
    % check the distance to c
    [d1km, d2km] = lldistkm(Xv(p, :), params.c);
    if d1km < params.R
        valid_idx(p) = 1; % add the sensor to valid list
        commMST(n_V+1, p) = d1km; % update comm. graph
    end
end

% start the greedy selection of best m_A sensors
for i = 1:params.m_A
    % print all valid indexes in this round
    fprintf('valid indexes in round %d:\n', i);
    for q = 1:length(valid_idx)
        if valid_idx(q) == 1
            fprintf('%d ', q);
        end
    end
    fprintf('\n');
    
    maxF = -Inf; % sensing quality of max result during searching
    maxM = -Inf; % maintenance cost of max result during searching
    maxRes = -Inf; % max result during searching
    maxRes_idx = -1; % the index of the best selection
    
    for j = 1:length(valid_idx)
        if valid_idx(j) == 1 && Xa_idx(j) == 0 % not select before
            % try to add this index to Xa
            Xa_idx(j) = 1;
            
            % calculate the current sensing quality
            X_remain = Xv(~Xa_idx, :);
            cov_remain = cov_vd(~Xa_idx, ~Xa_idx);
            Xa_cur = Xv(logical(Xa_idx), :);
            Ta_cur = Tv(logical(Xa_idx), :);
            cov_Xa_cur = cov_vd(logical(Xa_idx), logical(Xa_idx));
            
            % pass only the MST of selected sensors
            MST_idx = logical(vertcat(Xa_idx, [1]));
            commMST_cur = commMST(MST_idx, MST_idx);
            % combine all together
            curF = sense_quality(X_remain, cov_remain, Xa_cur, cov_Xa_cur, params.K);
            fprintf('sensing quality with node %d is %f\n', j, curF);
            curM = maintain_cost(Xa_cur, Ta_cur, commMST_cur, logging);
            fprintf('maintenance cost with node %d is %f\n', j, curM);
            curRes = params.alpha * (curF - lastF) + ...
                (1 - params.alpha) * (curM - lastM);
            
            % compare and update
            if curRes > maxRes && curM < params.Cm
                maxRes = curRes;
                maxRes_idx = j;
                maxF = curF;
                maxM = curM;
            end
            
            % reset the Xa index
            Xa_idx(j) = 0;
        end
    end
    
    % update the greedy selection
    if maxRes > -Inf
        Xa_idx(maxRes_idx) = 1; % add the sensors to the list
        lastF = maxF;
        lastM = maxM;
        fprintf('The selection in round %d is %d: [%f %f]\n', ...
            i, maxRes_idx, Xv(maxRes_idx, 1), Xv(maxRes_idx, 2));
        
        % update the valid indexes
        for k = 1:length(valid_idx)
            if valid_idx(k) == 1 % if is already valid, skip
                continue
            end
        
            % check the distance to the selected sensor in this round
            [d1km, d2km] = lldistkm(Xv(k, :), Xv(maxRes_idx, :));
            if d1km < params.R
                valid_idx(k) = 1; % add the sensor to valid list
                commMST(maxRes_idx, k) = d1km; % update comm. graph
            end
        end
    else
        error('No valid indexes to select!');
    end
end
% return the final selection solution
Xa = Xv(logical(Xa_idx), :);
% pass only the MST of selected sensors
MST_idx = logical(vertcat(Xa_idx, [1]));
commMST = commMST(MST_idx, MST_idx);
end

