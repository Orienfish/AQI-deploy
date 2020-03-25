function [Xa] = IDSQ(m_A, Xv, cov_vd, Tv, K, alpha, c, R)
%% Greedy heuristic IDSQ to place sensors on a subset of locations.
%
% Args:
%   m_A: number of sensors
%   Xv: a list of candidate locations to choose from
%   cov_vd: cov matrix at Xv given pre-deployment D
%   Tv: average temperature estimation at Xv in Celsius
%   K: the fitted RBF kernel function
%   alpha: the weight factor in IDSQ
%   c: position of the sink in [lat lon]
%   R: communication range of the sensors in km
%
% Return:
%   Xa: a list of solution locations to place sensors
addpath('./libs/');
addpath('./lldistkm/');
n_V = size(Xv, 1); % number of reference locations Xv
valid_idx = zeros(n_V, 1); % a list of valid indexes in comm. range
Xa_idx = zeros(n_V, 1); % init idx list for Xa to all zero
lastF = 0.0; % sensing quality in last round of greedy selection
commMST = NaN(n_V + 1); % init a matrix for connection graph
                        % first n entries are for Xv, the last one is for 
                        % the sink c. A single-direction MST.

% get the valid indexes directly connected to the sink
for p = 1:length(valid_idx)
    % check the distance to c
    [d1km, d2km] = lldistkm(Xv(p, :), c);
    if d1km < R
        valid_idx(p) = 1; % add the sensor to valid list
        commMST(n_V+1, p) = d1km; % update comm. graph
    end
end

% start the greedy selection of best m_A sensors
for i = 1:m_A
    % print all valid indexes in this round
    fprintf('valid indexes in round %d:\n', i);
    for q = 1:length(valid_idx)
        if valid_idx(q) == 1
            fprintf('%d ', q);
        end
    end
    fprintf('\n');
    
    maxF = 0.0; % max sensing quality during searching
    maxRes = 0.0; % max result during searching
    maxRes_idx = 0; % the index of the best selection
    
    for j = 1:length(valid_idx)
        if valid_idx(j) == 1 && Xa_idx(j) == 0 % not select before
            % try to add this index to Xa
            Xa_idx(j) = 1;
            
            % calculate the current sensing quality
            X_remain = Xv(~Xa_idx, :);
            fprintf('size of X_remain: %d\n', size(X_remain, 1));
            cov_remain = cov_vd(~Xa_idx, ~Xa_idx);
            fprintf('size of cov_remain: %d x %d\n', size(cov_remain, 1), ...
                size(cov_remain, 2));
            Xa_cur = Xv(logical(Xa_idx), :);
            Ta_cur = Tv(logical(Xa_idx), :);
            cov_Xa_cur = cov_vd(logical(Xa_idx), logical(Xa_idx));
            fprintf('size of X_cur: %d\n', size(Xa_cur, 1));
            fprintf('size of cov_Xa_cur: %d x %d\n', size(cov_Xa_cur, 1), ...
                size(cov_Xa_cur, 2));
            MST_idx = logical(vertcat(Xa_idx, [1]));
            commMST_cur = commMST(MST_idx, MST_idx);
            curF = sense_quality(X_remain, cov_remain, Xa_cur, cov_Xa_cur, K);
            fprintf('sensing quality gain at node %d is %f\n', j, curF - lastF);
            curRes = alpha * (curF - lastF) + (1 - alpha) * ...
                maintain_cost(Xa_cur, Ta_cur, commMST_cur);
            
            % compare and update
            if curRes > maxRes
                maxRes = curRes;
                maxRes_idx = j;
                maxF = curF;
            end
            
            % reset the Xa index
            Xa_idx(j) = 0;
        end
    end
    
    % update the greedy selection
    if maxRes > 0
        Xa_idx(maxRes_idx) = 1; % add the sensors to the list
        lastF = maxF;
        fprintf('The selection in round %d is %d: [%f %f]\n', ...
            i, maxRes_idx, Xv(maxRes_idx, 1), Xv(maxRes_idx, 2));
        
        % update the valid indexes
        for k = 1:length(valid_idx)
            if valid_idx(k) == 1 % if is already valid, skip
                continue
            end
        
            % check the distance to the selected sensor in this round
            [d1km, d2km] = lldistkm(Xv(k, :), Xv(maxRes_idx, :));
            if d1km < R
                valid_idx(k) = 1; % add the sensor to valid list
                commMST(maxRes_idx, k) = d1km; % update comm. graph
            end
        end
    else
        error('No valid indexes to select!');
    end
end
Xa = Xv(logical(Xa_idx), :);
end

