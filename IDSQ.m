function [res] = IDSQ(Qparams, params, IDSQparams)
%% Greedy heuristic IDSQ to place sensors on a subset of locations.
%
% Args:
%   Qparams.Xv: a list of candidate locations to choose from
%   Qparams.cov_vd: cov matrix at Qparams.Xv given pre-deployment D
%   Qparams.Xd: list of predeployment locations
%   Qparams.mean_temp_d: mean temperature at D
%   Qparams.cov_temp_d: cov matrix of temperature at D
%
%   params.m_A: number of sensors to deploy
%   params.Cm: maintenance cost budget
%   params.K: the fitted RBF kernel function
%   parmas.K_temp: the fitted RBF kernel function for temperature
%   params.c: position of the sink in [lat lon]
%   params.R: communication range of the sensors in km
%   params.logging: logging flag
%
%   IDSQparams.alpha: the weight factor in IDSQ
%
% Return:
%   res.Xa: a list of solution locations to place sensors
%   res.commMST: the generated communication graph
%   res.F: sensing quality of the solution
%   res.M: maintenance cost of the solution
addpath('./mlibs/');
addpath('./lldistkm/');
addpath('./gp/');

% local variables
n_V = size(Qparams.Xv, 1);  % number of reference locations Qparams.Xv
valid_idx = zeros(n_V, 1);  % a list of valid indexes in comm. range
Xa_idx = zeros(n_V, 1);     % init idx list for Xa to all zero
lastF = 0.0;                % sensing quality in last round of greedy selection
lastM = 0.0;                % maintenance cost in last round of greedy selection
commMST = NaN(n_V + 1);     % init a matrix for connection graph
                            % first n entries are for Qparams.Xv, the last one is  
                            % for the sink c. A single-direction MST.
predMST = NaN(n_V + 1, 1);  % predecessor nodes of the MST

% predict the ambient temperature at Xv in Celsius
[temp_mean_vd, temp_cov_vd] = gp_predict_knownD( ...
    Qparams.Xv, Qparams.Xd, Qparams.mean_temp_d, Qparams.cov_temp_d, ...
    params.K_temp);
temp_mean_vd = temp_mean_vd / 4 + 180; % weird fix
Tv = fah2cel(temp_mean_vd);  % convert to Celsius

% get the valid indexes directly connected to the sink
for p = 1:length(valid_idx)
    % check the distance to c
    [d1km, d2km] = lldistkm(Qparams.Xv(p, :), params.c);
    if d1km < params.R
        valid_idx(p) = 1; % add the sensor to valid list
        commMST(n_V+1, p) = d1km; % update comm. graph
        predMST(p) = n_V+1; % update the predecessor to sink
    end
end

% start the greedy selection of best m_A sensors
for i = 1:params.m_A
    % print all valid indexes in this round
    if params.logging
        fprintf('valid indexes in round %d:\n', i);
        for q = 1:length(valid_idx)
            if valid_idx(q) == 1
                fprintf('%d ', q);
            end
        end
        fprintf('\n');
    end
    
    maxF = -Inf;        % sensing quality of max result during searching
    maxM = -Inf;        % maintenance cost of max result during searching
    maxRes = -Inf;      % max result during searching
    maxRes_idx = -1;    % the index of the best selection
    
    for j = 1:length(valid_idx)
        if valid_idx(j) == 1 && Xa_idx(j) == 0 % not select before
            % try to add this index to Xa
            Xa_idx(j) = 1;
            
            % calculate the current sensing quality
            X_remain = Qparams.Xv(~Xa_idx, :);
            cov_remain = Qparams.cov_vd(~Xa_idx, ~Xa_idx);
            Xa_cur = Qparams.Xv(logical(Xa_idx), :);
            cov_Xa_cur = Qparams.cov_vd(logical(Xa_idx), logical(Xa_idx));
            
            % combine all together
            curF = sense_quality(X_remain, cov_remain, Xa_cur, ...
                cov_Xa_cur, params.K);
            curF = real(curF); % take the real part
            curM = maintain_cost(Qparams.Xv, Tv, Xa_idx, ...
                commMST, predMST, params.logging);
            curRes = IDSQparams.alpha * (curF - lastF) + ...
                (1 - IDSQparams.alpha) * (curM.C - lastM);
            
            if params.logging
                fprintf('node: %d senQ: %f main cost: %f res: %f\n', ...
                    j, curF, curM.C, curRes);
            end
            
            % compare and update
            if curRes > maxRes && curM.C < params.Cm
                maxRes = curRes;
                maxRes_idx = j;
                maxF = curF;
                maxM = curM.C;
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
        if params.logging
            fprintf('The selection in round %d is %d: [%f %f]\n', ...
                i, maxRes_idx, Qparams.Xv(maxRes_idx, 1), ...
                Qparams.Xv(maxRes_idx, 2));
        end
        
        % update the valid indexes
        for k = 1:length(valid_idx)
            if valid_idx(k) == 1 % if is already valid, skip
                continue
            end
        
            % check the distance to the selected sensor in this round
            [d1km, d2km] = lldistkm(Qparams.Xv(k, :), Qparams.Xv(maxRes_idx, :));
            if d1km < params.R
                valid_idx(k) = 1; % add the sensor to valid list
                commMST(maxRes_idx, k) = d1km; % update comm. graph
                predMST(k) = maxRes_idx; % update the predecessor
            end
        end
    else
        error('No valid indexes to select!');
    end
end
% return the final selection solution
res.Xa = Qparams.Xv(logical(Xa_idx), :);
% pass only the MST of selected sensors
MST_idx = logical(vertcat(Xa_idx, [1]));
res.commMST = commMST(MST_idx, MST_idx);
% pass the final results
res.F = lastF;
res.M = lastM;
end

