function [res] = DWG(Qparams, params)
%% Greedy heuristic IDSQ to place sensors on a subset of locations.
%
% Args:
%   Qparams.Xv: a list of candidate locations to choose from
%   Qparams.cov_vd: cov matrix at Qparams.Xv given pre-deployment D
%   Qparams.Xd: list of predeployment locations
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
%   params.logging: logging flag
%
% Return:
%   res.Xa: a list of solution locations to place sensors
%   res.commMST: the generated communication graph
%   res.F: sensing quality of the solution
%   res.M: maintenance cost of the solution
%   res.pred: the predecessors of each node

addpath('./mlibs/');
addpath('./lldistkm/');
addpath('./gp/');

% local variables
n_V = size(Qparams.Xv, 1);  % number of reference locations Qparams.Xv
valid_idx = zeros(n_V, 1);  % a list of valid indexes in comm. range
Xa_idx = zeros(n_V, 1);     % init idx list for Xa to all zero
lastF = 0.0;                % sensing quality in last round of greedy selection
commMST = NaN(n_V + 1);     % init a matrix for connection graph
                            % first n entries are for Qparams.Xv.
                            % A mutual-direction MST.
predMST = NaN(n_V + 1, 1);  % predecessor nodes of the MST

% get the valid indexes directly connected to the sink
predMST(n_V+1) = 0;         % configure the predecessor of the sink to 0
for p = 1:length(valid_idx)
    % check the distance to c
    [d1km, d2km] = lldistkm(Qparams.Xv(p, :), params.c);
    if d1km < params.R
        valid_idx(p) = 1;   % add the sensor to valid list
        commMST(n_V+1, p) = d1km; % update comm. graph
        predMST(p) = n_V+1; % update the predecessor to sink
    end
end

% create the undirected dist matrix for later use
D = zeros(params.n_V); % init a symmetric matrix for connection graph
nodes = Qparams.Xv;
for p = 1:params.n_V
    for q = p+1:params.n_V
        % check the distance to c
        [d1km, d2km] = lldistkm(nodes(p, :), nodes(q, :));
        % update the matrix
        D(p, q) = d1km;
        D(q, p) = d1km;
    end
end

%% create an information matrix based on all of the sensors, choose
% the two with biggest benefit-cost ratio as the base case
I = zeros(params.n_V);
for p = 1:params.n_V
    for q = p+1:params.n_V
        
        % calculate F(C1UC2)
        Xa_idx = zeros(n_V, 1);
        Xa_idx(p) = 1;
        Xa_idx(q) = 1;
        Xa_cur = Qparams.Xv(logical(Xa_idx), :);
        cov_Xa_cur = Qparams.cov_vd(logical(Xa_idx), logical(Xa_idx));
        fTotal = sense_quality(Qparams.Xv, Qparams.cov_vd, Xa_cur, ...
            cov_Xa_cur, params.K);
        
        % calcualte F(C1)
        Xa_idx(q) = 0;
        Xa_cur = Qparams.Xv(logical(Xa_idx), :);
        cov_Xa_cur = Qparams.cov_vd(logical(Xa_idx), logical(Xa_idx));
        f1 = sense_quality(Qparams.Xv, Qparams.cov_vd, Xa_cur, ...
            cov_Xa_cur, params.K);
        
        % calcualte F(C2)
        Xa_idx(p) = 0;
        Xa_idx(q) = 1;
        Xa_cur = Qparams.Xv(logical(Xa_idx), :);
        cov_Xa_cur = Qparams.cov_vd(logical(Xa_idx), logical(Xa_idx));
        f2 = sense_quality(Qparams.Xv, Qparams.cov_vd, Xa_cur, ...
            cov_Xa_cur, params.K);
        
        % min (fTotal - f(1 or 2)) => max f(i)
        fUse = f1;
        if f1 < f2
            fUse = f2;
        end
        
        numerator = fTotal - fUse; % get the numerator
        dist = D(p, q); % get the denominator
        I(p, q) = numerator / dist;
        I(q, p) = I(p, q);
    end
end

% get valid x and y coordinates for the biggest gain
M = 0.0; % largest element
K = [0, 0]; % index storing largest element
x = 0.0; % x-coordinate
y = 0.0; % y-coordinate

% continuously selecting biggest gain from the matrix; if not valid, set
% to -inf and select the next biggest gain from the matrix until valid
while true
    [M, K] = max(I(:)); % M: largest element; K: coordinate
    [x, y] = ind2sub(size(I),K); % get x and y coordinates of largest element
    if valid_idx(x) == 1 && valid_idx(y) == 1 % check valid
        break
    else
        I(x, y) = -inf; % not valid, set to -inf and start next round of selection
    end
end

% reset Xa_idx
Xa_idx = zeros(n_V, 1);
% add the sensors to the list
Xa_idx(x) = 1;
Xa_idx(y) = 1;
% getting current largest sensing quality
Xa_cur = Qparams.Xv(logical(Xa_idx), :);
cov_Xa_cur = Qparams.cov_vd(logical(Xa_idx), logical(Xa_idx));
lastF = sense_quality(Qparams.Xv, Qparams.cov_vd, Xa_cur, ...
            cov_Xa_cur, params.K);

%% start the greedy selection until satisfy the quota or m_A
cnt = 2;                    % the number of current selected sensors 
while true
    % print all valid indexes in this round
    if params.logging
        fprintf('valid indexes in round %d:\n', cnt);
        for q = 1:length(valid_idx)
            if valid_idx(q) == 1
                fprintf('%d ', q);
            end
        end
        fprintf('\n');
    end
    
    maxRes = -Inf;      % max result during searching
    maxRes_idx = -1;    % the index of the best selection
    % lastF: sensing quality in the last round
    
    for j = 1:n_V
        if valid_idx(j) == 1 && Xa_idx(j) == 0 % not select before
            % try to add this index to Xa
            Xa_idx(j) = 1;
            
            % calculate the current sensing quality (F(C1UC2))
            Xa_cur = Qparams.Xv(logical(Xa_idx), :);
            cov_Xa_cur = Qparams.cov_vd(logical(Xa_idx), logical(Xa_idx));
            curF = sense_quality(Qparams.Xv, Qparams.cov_vd, Xa_cur, ...
                cov_Xa_cur, params.K);
            curF = real(curF); % take the real part
            
            % calculate numerator
            % same logic as before: C1: selected cluster; C2: new node
            % F(C2)
            temp = zeros(n_V, 1);
            temp(j) = 1;
            temp_cur = Qparams.Xv(logical(temp), :);
            temp_cov = Qparams.cov_vd(logical(temp), logical(temp));
            nodeF = sense_quality(Qparams.Xv, Qparams.cov_vd, temp_cur, ...
                        temp_cov, params.K);
            nodeF = real(nodeF);
                    
            usedF = nodeF;
            if nodeF < lastF
                usedF = lastF;
            end
            
            numerator = curF - usedF;
            
            % calculate denominator
            % define the location of the cluster as the middle point
            % get the location of the current cluster
            totX = 0.0;
            totY = 0.0;
            for index = 1:n_V
                if Xa_idx(index) == 1 % selected points
                    if index ~= j
                        currentN = nodes(index, :);
                        totX = totX + currentN(1);
                        totY = totY + currentN(2);
                    end
                end
            end

            % get location of the cluster
            aveX = totX / cnt;
            aveY = totY / cnt;
            loc = [aveX aveY];
            % get distance
            [d1km, d2km] = lldistkm(loc, nodes(j, :));
            distance = d1km;
            
            % get gain
            curRes = numerator / distance;
            
            % compare and update
            if curRes > maxRes
                maxRes = curRes;
                maxRes_idx = j;
            end
            
            % reset the Xa index
            Xa_idx(j) = 0;
        end
    end
    
    % update the greedy selection
    if maxRes > -Inf
        Xa_idx(maxRes_idx) = 1; % add the sensors to the list
        lastF = maxRes;
        if params.logging
            fprintf('The selection in round %d is %d, senQ: %f\n', ...
                cnt, maxRes_idx, lastF);
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
    
    cnt = cnt + 1;
    
    % break the loop
    if lastF >= params.Q || cnt >= params.m_A
        break;
    end
end

% return the final selection solution
res.Xa = Qparams.Xv(logical(Xa_idx), :);
% pass only the MST of selected sensors
MST_idx = logical(vertcat(Xa_idx, [1]));
res.commMST = commMST(MST_idx, MST_idx);
% pass the final results
res.F = lastF;
% pass the predecossors
res.pred = predMST;

% calculate maintenance cost
% predict the ambient temperature at Xv in Celsius
[temp_mean_vd, temp_cov_vd] = gp_predict_knownD( ...
    Qparams.Xv, Qparams.Xd, Qparams.mean_temp_d, Qparams.cov_temp_d, ...
    params.K_temp);
Tv = fah2cel(temp_mean_vd); % convert to Celsius
res.M = maintain_cost(Qparams.Xv, Tv, Xa_idx, ...
                commMST, predMST, false);
end

