function [out] = ABC(Qparams, params, ABCparams)
%% Run the Artificial Bee Colony Optimization.
%
% Args:
%   Qparams.Xv: list of reference locations to predict, [lat lon]
%   Qparams.cov_vd: cov matrix at Xv given pre-deployment D
%   Qparams.Xd: list of predeployment locations
%   Qparams.mean_d: mean value at D
%   Qparams.cov_d: cov matrix at D
%   Qparams.Xa: list of locations we are supposed to observe, [lat lon]
%   Qparams.Ta: average temperature estimation at Xa in Celsius
%   Qparams.cov_ad: cov matrix at Xa given pre-deployment D
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
%   params.bound: bound for the area
%   params.PSOparams.weights: 1x3 vectors for the PSOparams.weights
%   params.penalty: penalty for non-connected nodes
%   params.logging: logging flag
%
%   ABCparams.nVar: number of unknown decision variables
%   ABCparams.VarSize: matrix size of decision variables
%   ABCparams.maxIter: maximum number of iterations
%   ABCparams.nPop: population size
%   ABCparams.nOnlooker: number of onlooker bees
%   ABCparams.L: abandonment limit parameter (trial limit)
%   ABCparams.a: acceleration coefficient upper bound
%   ABCparams.thres: penalty threshold in initialzation
%
% Return:
%   out.Position: the global best locations
%   out.Cost: the global best cost
%   out.senQuality: the global best sensing quality
%   out.mainCost: the global best maintenance cost
%   out.BestCosts: the global best cost iteration curve
addpath('./mlibs/');
addpath('./lldistkm/');
addpath('./gp/');

%% Initialization
% Empty Bee Structure
empty_bee.Position = [];
empty_bee.Cost = [];
% Initialize Population Array
pop = repmat(empty_bee, ABCparams.nPop, 1);
% Initialize Best Solution Ever Found
BestSol.Cost = inf;
% Note: no need to init Position, senQuality, mainCost
% they will be init in the initialization process

% Create Initial Population
for i = 1:ABCparams.nPop
    while true
        % generate random solution
        pop(i).Position(:, 1) = unifrnd(params.bound.latLower, ...
            params.bound.latUpper, ABCparams.nVar, 1);
        pop(i).Position(:, 2) = unifrnd(params.bound.lonLower, ...
            params.bound.lonUpper, ABCparams.nVar, 1);

        % cost evaluation
        [pm2_5_mean_ad, pm2_5_cov_ad] = gp_predict_knownD( ...
            pop(i).Position, Qparams.Xd, Qparams.mean_d, ...
            Qparams.cov_d, params.K);
        [temp_mean_ad, temp_cov_ad] = gp_predict_knownD( ...
            pop(i).Position, Qparams.Xd, Qparams.mean_temp_d, ...
            Qparams.cov_temp_d, params.K_temp);

        % setting the rest quality parameters
        Qparams.Xa = pop(i).Position;
        Qparams.Ta = fah2cel(temp_mean_ad);
        Qparams.cov_ad = pm2_5_cov_ad;
        res = costFunction(Qparams, params);
        
        % loop until generate one valid solutioin
        if res.P <= ABCparams.thres
            fprintf('init %d\n', i);
            break;
        end
    end
    
    pop(i).Cost = res.cost;
    pop(i).senQuality = res.F;
    pop(i).mainCost = res.M;
    
    % update the global best
    if pop(i).Cost <= BestSol.Cost
        BestSol = pop(i);
    end
end

% Abandonment Counter
C = zeros(ABCparams.nPop,1);

% Array to Hold Best Cost Values
BestCosts = zeros(ABCparams.maxIter, 1);


%% ABC Main Loop
for it = 1:ABCparams.maxIter
    % Recruited Bees
    for i = 1:ABCparams.nPop
        % Choose k randomly, not equal to i
        K = [1:i-1 i+1:ABCparams.nPop];
        k = K(randi([1 numel(K)]));
        
        % Define Acceleration Coeff.
        phi = ABCparams.a * unifrnd(-1, +1, ABCparams.VarSize);
        
        % New Bee Position
        newbee.Position = pop(i).Position + phi .* (pop(i).Position - ...
            pop(k).Position);
        %disp('position shift');
        %disp(phi .* (pop(i).Position - pop(k).Position));
        
        % apply lower and upper params.bound limits
        newbee.Position(:, 1) = ...
            max(newbee.Position(:, 1), params.bound.latLower);
        newbee.Position(:, 1) = ...
            min(newbee.Position(:, 1), params.bound.latUpper);
        newbee.Position(:, 2) = ...
            max(newbee.Position(:, 2), params.bound.lonLower);
        newbee.Position(:, 2) = ...
            min(newbee.Position(:, 2), params.bound.lonUpper);
        
        % Evaluation
        [pm2_5_mean_ad, pm2_5_cov_ad] = gp_predict_knownD( ...
            newbee.Position, Qparams.Xd, Qparams.mean_d, ...
            Qparams.cov_d, params.K);
        [temp_mean_ad, temp_cov_ad] = gp_predict_knownD( ...
            newbee.Position, Qparams.Xd, Qparams.mean_temp_d, ...
            Qparams.cov_temp_d, params.K_temp);

        % setting the rest quality parameters
        Qparams.Xa = newbee.Position;
        Qparams.Ta = fah2cel(temp_mean_ad);
        Qparams.cov_ad = pm2_5_cov_ad;
        res = costFunction(Qparams, params);
        
        newbee.Cost = res.cost;
        newbee.senQuality = res.F;
        newbee.mainCost = res.M;
        
        % Comparision
        if newbee.Cost <= pop(i).Cost
            pop(i) = newbee;
        else
            C(i) = C(i) + 1;
        end     
    end
    
    % Calculate Fitness Values and Selection Probabilities
    F = zeros(ABCparams.nPop, 1);
    MeanCost = mean([pop.Cost]);
    for i = 1:ABCparams.nPop
        F(i) = exp(-pop(i).Cost / MeanCost); % Convert Cost to Fitness
    end
    P = F / sum(F);
    
    % Onlooker Bees
    for m=1:ABCparams.nOnlooker
        % Select Source Site
        i = RouletteWheelSelection(P);
        
        % Choose k randomly, not equal to i
        K = [1:i-1 i+1:ABCparams.nPop];
        k = K(randi([1 numel(K)]));
        
        % Define Acceleration Coeff.
        phi = ABCparams.a * unifrnd(-1, +1, ABCparams.VarSize);
        
        % New Bee Position
        newbee.Position = pop(i).Position + phi .* (pop(i).Position - ...
            pop(k).Position);
        
        % apply lower and upper params.bound limits
        newbee.Position(:, 1) = ...
            max(newbee.Position(:, 1), params.bound.latLower);
        newbee.Position(:, 1) = ...
            min(newbee.Position(:, 1), params.bound.latUpper);
        newbee.Position(:, 2) = ...
            max(newbee.Position(:, 2), params.bound.lonLower);
        newbee.Position(:, 2) = ...
            min(newbee.Position(:, 2), params.bound.lonUpper);
        
        % Evaluation
        [pm2_5_mean_ad, pm2_5_cov_ad] = gp_predict_knownD( ...
            newbee.Position, Qparams.Xd, Qparams.mean_d, ...
            Qparams.cov_d, params.K);
        [temp_mean_ad, temp_cov_ad] = gp_predict_knownD( ...
            newbee.Position, Qparams.Xd, Qparams.mean_temp_d, ...
            Qparams.cov_temp_d, params.K_temp);

        % setting the rest quality parameters
        Qparams.Xa = newbee.Position;
        Qparams.Ta = fah2cel(temp_mean_ad);
        Qparams.cov_ad = pm2_5_cov_ad;
        res = costFunction(Qparams, params);
        
        newbee.Cost = res.cost;
        newbee.senQuality = res.F;
        newbee.mainCost = res.M;
        
        % Comparision
        if newbee.Cost <= pop(i).Cost
            pop(i) = newbee;
        else
            C(i) = C(i)+1;
        end
    end
    
    % Scout Bees
    for i = 1:ABCparams.nPop
        if C(i) >= ABCparams.L
            fprintf('Abandon current source!\n');
            % generate random solution
            pop(i).Position(:, 1) = unifrnd(params.bound.latLower, ...
                params.bound.latUpper, ABCparams.nVar, 1);
            pop(i).Position(:, 2) = unifrnd(params.bound.lonLower, ...
                params.bound.lonUpper, ABCparams.nVar, 1);
            
            % cost evaluation
            [pm2_5_mean_ad, pm2_5_cov_ad] = gp_predict_knownD( ...
                pop(i).Position, Qparams.Xd, Qparams.mean_d, ...
                Qparams.cov_d, params.K);
            [temp_mean_ad, temp_cov_ad] = gp_predict_knownD( ...
                pop(i).Position, Qparams.Xd, Qparams.mean_temp_d, ...
                Qparams.cov_temp_d, params.K_temp);

            % setting the rest quality parameters
            Qparams.Xa = pop(i).Position;
            Qparams.Ta = fah2cel(temp_mean_ad);
            Qparams.cov_ad = pm2_5_cov_ad;
            res = costFunction(Qparams, params);
            
            pop(i).Cost = res.cost;
            pop(i).senQuality = res.F;
            pop(i).mainCost = res.M;
            
            % reset counter
            C(i) = 0;
        end
    end
    
    % Update Best Solution Ever Found
    for i = 1:ABCparams.nPop
        if pop(i).Cost < BestSol.Cost
            %fprintf('update best solution!\n');
            BestSol = pop(i);
        end
    end
    
    % Store Best Cost Ever Found
    BestCosts(it)=BestSol.Cost;
    
    % Display Iteration Information
    fprintf('Iteration %d: Best Cost: %f senQ: %f mainCost: %f\n', ...
        it, BestSol.Cost, BestSol.senQuality, BestSol.mainCost);
end

out = BestSol; % include position, cost, senQuality, mainCost
out.BestCosts = BestCosts; % cost iteration curve
end

