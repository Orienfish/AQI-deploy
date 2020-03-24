function [Xa] = IDSQ(m_A, Xv, cov_vd, K, alpha)
%% Greedy heuristic IDSQ to place sensors on a subset of locations.
%
% Args:
%   m_A: number of sensors
%   Xv: a list of candidate locations to choose from
%   cov_vd: cov matrix at Xv given pre-deployment D
%   K: the fitted RBF kernel function
%   alpha: the weight factor in IDSQ
%
% Return:
%   Xa: a list of solution locations to place sensors
addpath('./libs/');
addpath('./lldistkm/');
valid_idx = zeros(1, size(Xv, 1)); % a list of valid indexes in comm. range
Xa_idx = zeros(1, size(Xv, 1)); % init idx list for Xa to all zero
lastF = 0.0; % sensing quality in last round of greedy selection

% start the greedy selection of best m_A sensors
for i = 1:m_A
    % get the valid indexes
    for j = 1:length(valid_idx)
        % check the distance to c
        [d1km, d2km] = lldistkm(Xv(j, :), c);
        if d1km < R
            valid_idx(j) = 1;
            continue;
        end
        
        % check the distance to selected sensors
        for k = 1:length(Xa_idx)
            if Xa_idx(k) == 1 % selected sensros
                [d1km, d2km] = lldistkm(Xv(j, :), Xv(k, :));
                if d1km < R:
                    valid_idx(k) = 1;
                    continue;
                end
            end
        end
    end
    % print all valid indexes in this round
    fprintf('valid indexes in round %d:\n', i);
    for j = 1:length(valid_idx)
        if valid_idx(j) == 1
            fprintf('%d ', j);
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
            cov_Xa = cov_vd(Xa_idx, Xa_idx);
            fprintf('size of cov_Xa: %d x %d\n', size(cov_Xa, 1), ...
                size(cov_Xa, 2));
            curF = sense_Quality(X_candidate, cov_candidate, Xa, cov_Xa, K);
            curRes = alpha * (curF - lastF) + (1 - alpha) * maintainCost(Xa);
            
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
    else
        error('No valid indexes to select!');
    end
end
Xa = Xv(Xa_idx, :);
end

