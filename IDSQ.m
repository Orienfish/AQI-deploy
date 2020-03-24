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
valid_idx = zeros(1, n_V); % a list of valid indexes in comm. range
Xa_idx = zeros(1, n_V); % init idx list for Xa to all zero
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

% start the greedy selection of best m_A sensors
for i = 1:m_A
    % print all valid indexes in this round
    fprintf('valid indexes in round %d:\n', i);
    for a = 1:length(valid_idx)
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
            X_candidate = Xv(~Xa_idx, :);
            fprintf('size of X_candidate: %d\n', size(X_candidate, 1));
            cov_candidate = cov_vd(~Xa_idx, ~Xa_idx);
            fprintf('size of cov_candidate: %d x %d\n', size(cov_candidate, 1), ...
                size(cov_candidate, 2));
            Xa = Xv(Xa_idx, :);
            Ta = Tv(Xa_idx, :);
            cov_Xa = cov_vd(Xa_idx, Xa_idx);
            fprintf('size of cov_Xa: %d x %d\n', size(cov_Xa, 1), ...
                size(cov_Xa, 2));
            curF = sense_Quality(X_candidate, cov_candidate, Xa, cov_Xa, K);
            curRes = alpha * (curF - lastF) + (1 - alpha) * maintain_cost(Xa, Ta, commMST);
            
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
        fprintf('The selection in round %d is %d\n', i, maxRes_idx);
        
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
Xa = Xv(Xa_idx, :);
end

