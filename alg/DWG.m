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
valid_idx = zeros(n_V, n_V);  % a matrix of valid indexes in comm. range
Xa_idx = zeros(n_V, 1);     % init idx list for Xa to all zero
clusters = zeros(n_V, n_V); % matrix storing information of clusters
                            % initially n_V by n_V matrix and each start
                            % with one cluster, gradually truncating the
                            % size based on merging
dist = zeros(params.n_V);   % init a symmetric matrix for connection graph
                            % build a distance matrix based on distance 
                            % between each nodes
nodes = Qparams.Xv;         % short-cut for the full list of nodes

% initialize the clusters matrix
for p = 1:n_V
    clusters(p, 1) = p;
end


% traverse through all the nodes and build the distance matrix and the
% valid_idx matrix
for p = 1:params.n_V
    for q = p+1:params.n_V
        % check the distance to c
        [d1km, d2km] = lldistkm(nodes(p, :), nodes(q, :));
        % update the dist matrix
        dist(p, q) = d1km;
        dist(q, p) = d1km;
        % update valid_idx matrix
        % if the distance between two nodes is smaller than communication
        % range, then it's valid to combine the two nodes
        if d1km < params.R
            valid_idx(p, q) = 1;
            valid_idx(q, p) = 1;
        end
    end
end

% dist: distance matrix between each pair of nodes
% D: distance matrix between each pair of clusters
% originally each cluster is one node, so D == dist
% this matrix will be modified (truncated) in the merging process
D = dist;

%% build sensing_quality matrix for later use
% this matrix will be modified (truncated) in the merging process
sq = zeros(params.n_V, params.n_V);
for p = 1:params.n_V
    % q not start from (p+1) bc need to consider the on-cluster case when
    % calculating gain
    for q = p:params.n_V
        % add all nodes in those two clusters into Xa_idx
        Xa_idx = add_to_Xa(clusters, p, q, n_V);
        % calculate the current sensing quality
        curF = cal_sq(Xa_idx, Qparams.Xv, Qparams.cov_vd, params.K);

        % update the sensing_quality matrix
        sq(p, q) = curF;
        sq(q, p) = curF;
    end
end

%% build information matrix storing gain if combining the two clusters
% this matrix will be modified (truncated) in the merging process
I = NaN(params.n_V, params.n_V);
for p = 1:params.n_V
    for q = p+1:params.n_V
        % skip if not valid
        if valid_idx(p, q) == 0
            continue
        end
        
        % get the gain if combining p and q
        g = cal_gain(sq, p, q, D);

        % update the matrix
        I(p, q) = g;
        I(q, p) = g;
    end
end
    
% directly begin merging process
while true
    numOfC = size(clusters); % number of clusters
    numOfC = numOfC(1); % number of rows is the number of clusters
    
    % check for number of clusters
    % exit the loop if combining all nodes
    if numOfC == 1
        Xa_idx = ones(n_V, 1);
        break
    end
    
%     needbreak = 0;
%     % check if m_A nodes have been selected
%     for i = 1:numOfC
%         if clusters(i, params.m_A) == 0
%             continue
%         end
%         fprintf("here");
%         % exit the loop
%         Xa_idx = zeros(n_V, 1);
%         for p = 1:n_V
%             if clusters(i, p) == 0
%                 break
%             end
%             Xa_idx(clusters(i, p)) = 1;
%         end
%         needbreak = 1;
%         break
%     end
%     
%     if needbreak == 1
%         break
%     end

    %% begin selection process
    % get valid x and y coordinates for the biggest gain (from matrix I,
    % information)
    ret = max_val(I, valid_idx);
    x = ret.x;
    y = ret.y;
    
    % find the first index where cluster(x, idx) = 0 for the combination
    % process: add all nodes in y to x
    fstIdx = 0;
    for i = 1:n_V
        if clusters(x, i) == 0
            fstIdx = i;
            break
        end
    end
    
    % put all nodes in y into x
    for i = 1:n_V
        % traversing all elements in y
        if clusters(y, i) == 0
            break
        end
        node = clusters(y, i);
        clusters(x, fstIdx) = node;
        fstIdx = fstIdx + 1; % update the index
    end
    
    %% update D, valid_idx, sq, and I matrices accordingly
    % D: matrix storing distance between clusters
    % sq: matrix storing sensing quality of the two clusters
    % I: matrix storing information (gain) of the two clusters
    
    % update D and valid_idx
    % update xth column and xth row
    for p = 1:numOfC
        if p == y
            continue
        end
        % shortest distance between two clusters
        minD = Inf;
        for i = 1:n_V
            % if reaching zero, it's the end of this cluster
            if clusters(x, i) == 0
                break
            end
            for j = 1:n_V
                % same logic: reaching 0 == reaching the end of this
                % cluster
                if clusters(p, j) == 0
                    break
                end
                % update minD if needed
                n1 = clusters(x, i);
                n2 = clusters(p, j);
                if dist(n1, n2) < minD
                    minD = dist(n1, n2);
                end
            end
        end
        % after finding minD, update
        D(x, p) = minD;
        D(p, x) = minD;
        
        % update valid_idx
        if minD < params.R
            valid_idx(p, x) = 1;
            valid_idx(x, p) = 1;
        end
    end
    
    % update sensing quality
    % update xth column and xth row
    for p = 1:numOfC
        if p == y
            continue
        end
        % add all nodes in those two clusters into Xa_idx
        Xa_idx = add_to_Xa(clusters, p, x, n_V);
        % calculate the current sensing quality
        curF = cal_sq(Xa_idx, Qparams.Xv, Qparams.cov_vd, params.K);

        sq(p, x) = curF;
        sq(x, p) = curF;
    end
    
    % update information matrix
    % update xth column and xth row
    for p = 1:numOfC
        if p == y
            continue
        end
        % skip if not valid
        if valid_idx(p, x) == 0
            continue
        end
        
        % get the gain if combining p and x
        g = cal_gain(sq, p, x, D);
        
        I(p, x) = g;
        I(x, p) = g;
    end
    
    % delete yth row in clusters
    clusters(y, :) = [];
    % delete yth column and yth row in valid_idx
    valid_idx(y, :) = [];
    valid_idx(:, y) = [];
    % delete yth column and yth row in D (distance matrix of clusters)
    D(y, :) = [];
    D(:, y) = [];
    % delete yth column and yth row in I (gain matrix)
    I(y, :) = [];
    I(:, y) = [];
    % delete yth column and yth row in sq (sensing quality matrix)
    sq(y, :) = [];
    sq(:, y) = [];
    
    %% check whether to end loop
    % biggest sensing quality
    ret = max_val(sq, valid_idx);
    M = ret.M;
    
    % if reaching the quota, exit the loop
    if M >= params.Q
        % update selectd nodes
        x = ret.x;
        y = ret.y;
        Xa_idx = add_to_Xa(clusters, x, y, n_V);
        break
    end
end


%% return values
res.Xa = Qparams.Xv(logical(Xa_idx), :);
% pass only the MST of selected sensors
[G, predMST] = MST(res.Xa, params.c, params.R);
res.commMST = G;
commMST = G;

% pass the predecossors
res.pred = predMST;

% calculate sensing quality
res.F = cal_sq(Xa_idx, Qparams.Xv, Qparams.cov_vd, params.K);

% calculate maintenance cost
% predict the ambient temperature at Xv in Celsius
[temp_mean_vd, temp_cov_vd] = gp_predict_knownD( ...
    Qparams.Xv, Qparams.Xd, Qparams.mean_temp_d, Qparams.cov_temp_d, ...
    params.K_temp);
Tv = fah2cel(temp_mean_vd); % convert to Celsius

res.M = maintain_cost(res.Xa, Tv, ~isnan(predMST), ...
                commMST, predMST, false);
end

%% the function to calculate sensing quality based on Xa_idx
function [curF] = cal_sq(Xa_idx, Xv, cov_vd, K)
    X_remain = Xv(~Xa_idx, :);
    cov_remain = cov_vd(~Xa_idx, ~Xa_idx);
    Xa_cur = Xv(logical(Xa_idx), :);
    cov_Xa_cur = cov_vd(logical(Xa_idx), logical(Xa_idx));
    curF = sense_quality(X_remain, cov_remain, Xa_cur, cov_Xa_cur, K);
    curF = real(curF); % take the real part
end

%% the function finds the max element in the matrix and also returns the index
function [ret] = max_val(matrix, valid_idx)
    % x: x coordinate of max element
    % y: y coordinate of max element
    % M: largest element
    [M, K] = max(matrix(logical(valid_idx))); % M: largest element; K: coordinate
    K = find(matrix == M);
    K = K(1);
    [x, y] = ind2sub(size(matrix), K); % get x and y coordinates of largest element
    ret.x = x;
    ret.y = y;
    ret.M = M;
end

%% the function adds nodes to Xa_idx and returns the new Xa_idx
function [Xa] = add_to_Xa(clusters, p, q, n_V)
    Xa_idx = zeros(n_V, 1);
    for i = 1:n_V
        % add to Xa_idx if not reaching 0
        if clusters(p, i) ~= 0
            currentN = clusters(p, i);
            Xa_idx(currentN) = 1;
        end
        if clusters(q, i) ~= 0
            currentN = clusters(q, i);
            Xa_idx(currentN) = 1;
        end
    end
    Xa = Xa_idx;
end

%% the function calculates gain
function [gain] = cal_gain(sq, p, q, D)
    % calculate numerator: F(C1 U C2) - max(F(Ci))
    % get fTotal: F(C1 U C2)
    fTotal = sq(p, q);

    % get F(C1) and F(C2)
    f1 = sq(p, p);
    f2 = sq(q, q);
    % get second part: max F(Ci)
    sndP = f1;
    if f2 > f1
        sndP = f2;
    end

    % get numerator
    numerator = fTotal - sndP;

    % get denominator, the distance between those two clusters
    denominator = D(p, q);
    
    gain = numerator / denominator;
end
