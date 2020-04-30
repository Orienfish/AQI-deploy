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
commMST = NaN(n_V + 1);     % init a matrix for connection graph
                            % first n entries are for Qparams.Xv.
                            % A mutual-direction MST.
predMST = NaN(n_V + 1, 1);  % predecessor nodes of the MST
clusters = zeros(n_V, n_V);   % matrix storing information of clusters
                            % initially n_V by n_V matrix and each start
                            % with one cluster, gradually truncating the
                            % size based on merging
for p = 1:n_V
    clusters(p, 1) = p;
end
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

% build a distance matrix based on distance between each nodes
% create the undirected dist matrix for later use
dist = zeros(params.n_V); % init a symmetric matrix for connection graph
nodes = Qparams.Xv;
for p = 1:params.n_V
    for q = p+1:params.n_V
        % check the distance to c
        [d1km, d2km] = lldistkm(nodes(p, :), nodes(q, :));
        % update the matrix
        dist(p, q) = d1km;
        dist(q, p) = d1km;
    end
end

% build distance matrix, each element is the distance between two clusters
% this matrix will be modified (truncated) in the merging process
D = dist;

%% build valid_idx matrix storing whether it's valid to combine the two clusters
% this matrix will be modified (truncated) in the merging process
valid_idx = zeros(params.n_V, params.n_V);
for p = 1:params.n_V
    for q = p+1:params.n_V
        if D(p, q) < params.R
            valid_idx(p, q) = 1;
            valid_idx(q, p) = 1;
        end
    end
end

%% build sensing_quality matrix for later use
% this matrix will be modified (truncated) in the merging process
sq = zeros(params.n_V, params.n_V);
for p = 1:params.n_V
    % q not start from (p+1) bc need to consider one cluster
    for q = p:params.n_V
        % add all nodes in those two clusters into Xa_idx
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
        % calculate the current sensing quality
        X_remain = Qparams.Xv(~Xa_idx, :);
        cov_remain = Qparams.cov_vd(~Xa_idx, ~Xa_idx);
        Xa_cur = Qparams.Xv(logical(Xa_idx), :);
        cov_Xa_cur = Qparams.cov_vd(logical(Xa_idx), logical(Xa_idx));
        curF = sense_quality(X_remain, cov_remain, Xa_cur, ...
            cov_Xa_cur, params.K);
        curF = real(curF); % take the real part

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
        % calculate numerator
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

        % get denominator
        denominator = D(p, q);

        I(p, q) = numerator / denominator;
        I(q, p) = I(p, q);
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

    %% begin selection process
    % get valid x and y coordinates for the biggest gain
    M = 0.0; % largest element
    K = [0, 0]; % index storing largest element
    x = 0.0; % xth cluster
    y = 0.0; % yth cluster
    [M, K] = max(I(logical(valid_idx))); % M: largest element; K: coordinate
    K = find(I == M);
    K = K(1);
    [x, y] = ind2sub(size(I),K); % get x and y coordinates of largest element
    
    % find the first index where cluster(x, idx) = 0
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
    % update D
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
    end
    
    % update valid_idx
    % update xth column and xth row
    for p = 1:numOfC
        if p == y
            continue
        end
        if D(p, x) < params.R
            valid_idx(p, x) = 1;
            valid_idx(x, p) = 1;
        end
    end
    
    % update sq
    % update xth column and xth row
    for p = 1:numOfC
        if p == y
            continue
        end
        % add all nodes in those two clusters into Xa_idx
        Xa_idx = zeros(n_V, 1);
        for i = 1:n_V
            % add to Xa_idx if not reaching 0
            if clusters(p, i) ~= 0
                currentN = clusters(p, i);
                Xa_idx(currentN) = 1;
            end
            if clusters(x, i) ~= 0
                currentN = clusters(x, i);
                Xa_idx(currentN) = 1;
            end
        end
        % calculate the current sensing quality
        X_remain = Qparams.Xv(~Xa_idx, :);
        cov_remain = Qparams.cov_vd(~Xa_idx, ~Xa_idx);
        Xa_cur = Qparams.Xv(logical(Xa_idx), :);
        cov_Xa_cur = Qparams.cov_vd(logical(Xa_idx), logical(Xa_idx));
        curF = sense_quality(X_remain, cov_remain, Xa_cur, ...
            cov_Xa_cur, params.K);
        curF = real(curF); % take the real part

        sq(p, x) = curF;
        sq(x, p) = curF;
    end
    
    % update I
    % update xth column and xth row
    for p = 1:numOfC
        if p == y
            continue
        end
        % skip if not valid
        if valid_idx(p, x) == 0
            continue
        end
        % calculate numerator
        % get fTotal: F(C1 U C2)
        fTotal = sq(p, x);
        
        % get F(C1) and F(C2)
        f1 = sq(p, p);
        f2 = sq(x, x);
        % get second part: max F(Ci)
        sndP = f1;
        if f2 > f1
            sndP = f2;
        end

        % get numerator
        numerator = fTotal - sndP;

        % get denominator
        denominator = D(p, x);
        
        I(p, x) = numerator / denominator;
        I(x, p) = I(p, x);
    end
    
    % delete yth row in clusters
    clusters(y, :) = [];
    % delete yth column and yth row in valid_idx
    valid_idx(y, :) = [];
    valid_idx(:, y) = [];
    % delete yth column and yth row in D
    D(y, :) = [];
    D(:, y) = [];
    % delete yth column and yth row in I
    I(y, :) = [];
    I(:, y) = [];
    % delete yth column and yth row in sq
    sq(y, :) = [];
    sq(:, y) = [];
    
    numOfC = numOfC - 1;
    
    %% check whether to end loop
    % biggest sensing quality
    M = 0.0; % largest element
    K = [0, 0]; % index storing largest element
    x = 0.0; % x-coordinate
    y = 0.0; % y-coordinate
    [M, K] = max(sq(logical(valid_idx))); % M: largest element; K: coordinate
    % if reaching the quota, exit the loop
    if M >= params.Q
        % update selectd nodes
        Xa_idx = zeros(n_V, 1);
        K = find(sq == M);
        K = K(1);
        [x, y] = ind2sub(size(sq),K);
        for i = 1:n_V
            if clusters(x, i) ~= 0
                currentN = clusters(x, i);
                Xa_idx(currentN) = 1;
            end
            if clusters(y, i) ~= 0
                currentN = clusters(y, i);
                Xa_idx(currentN) = 1;
            end
        end
        break
    end
end


%% return values
res.Xa = Qparams.Xv(logical(Xa_idx), :);
% pass only the MST of selected sensors
MST_idx = logical(vertcat(Xa_idx, [1]));
res.commMST = commMST(MST_idx, MST_idx);

% pass the predecossors
res.pred = predMST;

% calculate sensing quality
X_remain = Qparams.Xv(~Xa_idx, :);
cov_remain = Qparams.cov_vd(~Xa_idx, ~Xa_idx);
Xa_cur = Qparams.Xv(logical(Xa_idx), :);
cov_Xa_cur = Qparams.cov_vd(logical(Xa_idx), logical(Xa_idx));
curF = sense_quality(X_remain, cov_remain, Xa_cur, ...
            cov_Xa_cur, params.K);
res.F = real(curF);

% calculate maintenance cost
% predict the ambient temperature at Xv in Celsius
[temp_mean_vd, temp_cov_vd] = gp_predict_knownD( ...
    Qparams.Xv, Qparams.Xd, Qparams.mean_temp_d, Qparams.cov_temp_d, ...
    params.K_temp);
Tv = fah2cel(temp_mean_vd); % convert to Celsius

%res.M = maintain_cost(Qparams.Xv, Tv, Xa_idx, ...
%                commMST, predMST, false);
res.M = 1;
end


