function [out] = PSO(Qparams, params, PSOparams)
%% Run the Particle SPSOparams.warm Optimization.
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
%   params.m_A: number of sensors to deploy
%   params.Cm: maintenance cost budget
%   params.K: the fitted RBF kernel function
%   parmas.K_temp: the fitted RBF kernel function for temperature
%   params.c: position of the sink in [lat lon]
%   params.R: communication range of the sensors in km
%   params.PSOparams.weights: 1x3 vectors for the PSOparams.weights
%   params.penalty: penalty for non-connected nodes
%   params.logging: logging flag
%
%   PSOparams.nVar: number of unknoPSOparams.wn decision variables
%   PSOparams.VarSize: matrix size of decision variables
%   PSOparams.maxIter: maximum number of iterations
%   PSOparams.nPop: population size
%   PSOparams.chi: constriction factor
%   PSOparams.w: inertia coefficient
%   PSOparams.wdamp: damping ratio of inertia coefficient
%   PSOparams.c1: personal acceleration coefficient
%   PSOparams.c2: social acceleration coefficient
%   PSOparams.bound: PSOparams.bound for the area
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

%% parameters of PSO
% velocity limits
MaxVelocity = 0.2 * [(PSOparams.bound.latUpper - PSOparams.bound.latLower), ...
    (PSOparams.bound.lonUpper - PSOparams.bound.lonLower)];
MinVelocity = - MaxVelocity;

%% initialization
% the particle template
empty_particle.Position = [];
empty_particle.Velocity = [];
empty_particle.Cost = [];
empty_particle.senQuality = [];
empty_particle.mainCost = [];
empty_particle.Best.Position = []; % personal best
empty_particle.Best.Cost = [];     % personal best
empty_particle.Best.senQuality = [];
empty_particle.Best.mainCost = [];

% create population array
particle = repmat(empty_particle, PSOparams.nPop, 1);

% initialize global best
GlobalBest.Cost = inf;
% Note: no need to init Position, senQuality, mainCost
% they PSOparams.will be init in the initialization process

% initialize population members
for i=1:PSOparams.nPop
    % generate random solution
    particle(i).Position(:, 1) = unifrnd(PSOparams.bound.latLower, ...
        PSOparams.bound.latUpper, PSOparams.nVar, 1);
    particle(i).Position(:, 2) = unifrnd(PSOparams.bound.lonLower, ...
        PSOparams.bound.lonUpper, PSOparams.nVar, 1);
    
    % initialize velocity
    particle(i).Velocity = zeros(PSOparams.VarSize);
    
    % cost evaluation
    [pm2_5_mean_ad, pm2_5_cov_ad] = gp_predict_knownD( ...
        particle(i).Position, Qparams.Xd, Qparams.mean_d, ...
        Qparams.cov_d, params.K);
    [temp_mean_ad, temp_cov_ad] = gp_predict_knownD( ...
        particle(i).Position, Qparams.Xd, Qparams.mean_temp_d, ...
        Qparams.cov_temp_d, params.K_temp);
    temp_mean_ad = temp_mean_ad / 4 + 180; % weird fix
    
    % setting the rest quality parameters
    Qparams.Xa = particle(i).Position;
    Qparams.Ta = fah2cel(temp_mean_ad);
    Qparams.cov_ad = pm2_5_cov_ad;
    res = costFunction(Qparams, params);
    
    particle(i).Cost = res.cost;
    particle(i).senQuality = res.F;
    particle(i).mainCost = res.M;
    
    % update the personal best
    particle(i).Best.Position = particle(i).Position;
    particle(i).Best.Cost = particle(i).Cost;
    particle(i).Best.senQuality = particle(i).senQuality;
    particle(i).Best.mainCost = particle(i).mainCost;
    
    % update the global best
    if particle(i).Best.Cost < GlobalBest.Cost
        GlobalBest = particle(i).Best;
    end
end

% array to hold best cost value on each iteration
BestCosts = zeros(PSOparams.maxIter, 1);

%% main loop of PSO
for it = 1:PSOparams.maxIter
    for i = 1:PSOparams.nPop
        % update velocity
        particle(i).Velocity = particle(i).Velocity ...
            + PSOparams.c1*rand(PSOparams.VarSize).*(particle(i).Best.Position - particle(i).Position) ...
            + PSOparams.c2*rand(PSOparams.VarSize).*(GlobalBest.Position - particle(i).Position);
        
        % apply velocity limits
        particle(i).Velocity = max(particle(i).Velocity, MinVelocity);
        particle(i).Velocity = min(particle(i).Velocity, MaxVelocity);
        
        % update position
        particle(i).Position = particle(i).Position + particle(i).Velocity;
        
        % apply loPSOparams.wer and upper PSOparams.bound limits
        particle(i).Position(:, 1) = ...
            max(particle(i).Position(:, 1), PSOparams.bound.latLower);
        particle(i).Position(:, 1) = ...
            min(particle(i).Position(:, 1), PSOparams.bound.latUpper);
        particle(i).Position(:, 2) = ...
            max(particle(i).Position(:, 2), PSOparams.bound.lonLower);
        particle(i).Position(:, 2) = ...
            min(particle(i).Position(:, 2), PSOparams.bound.lonUpper);
        
        % cost evaluation
        [pm2_5_mean_ad, pm2_5_cov_ad] = gp_predict_knownD( ...
            particle(i).Position, Qparams.Xd, Qparams.mean_d, ...
            Qparams.cov_d, params.K);
        [temp_mean_ad, temp_cov_ad] = gp_predict_knownD( ...
            particle(i).Position, Qparams.Xd, Qparams.mean_temp_d, ...
            Qparams.cov_temp_d, params.K_temp);
        temp_mean_ad = temp_mean_ad / 4 + 180; % weird fix
        
        % Qparams.Xv = V;                 % use the same value as init
        % Qparams.cov_vd = pm2_5_cov_vd;  % use the same value as init
        Qparams.Xa = particle(i).Position;
        Qparams.Ta = fah2cel(temp_mean_ad);
        Qparams.cov_ad = pm2_5_cov_ad;
        res = costFunction(Qparams, params);

        particle(i).Cost = res.cost;
        particle(i).senQuality = res.F;
        particle(i).mainCost = res.M;
        
        if params.logging
            fprintf('particle %d\n', i);
            fprintf('Cost: %f senQ: %f mainCost: %f\n', particle(i).Cost, ...
                particle(i).senQuality, particle(i).mainCost);
        end
        
        % update personal best
        if particle(i).Cost < particle(i).Best.Cost
            fprintf('update personal best\n');
            particle(i).Best.Position = particle(i).Position;
            particle(i).Best.Cost = particle(i).Cost;
            particle(i).Best.senQuality = particle(i).senQuality;
            particle(i).Best.mainCost = particle(i).mainCost;
            
            % update the global best
            if particle(i).Best.Cost < GlobalBest.Cost
                fprintf('update global best\n');
                GlobalBest = particle(i).Best;
            end
        end
    end
    
    % store the best cost value
    BestCosts(it) = GlobalBest.Cost;
    
    % display iteration information
    fprintf('Iteration %d: Best Cost: %f senQ: %f mainCost: %f\n', ...
        it, GlobalBest.Cost, GlobalBest.senQuality, GlobalBest.mainCost);
    
    % damping inertia coefficient
    PSOparams.w = PSOparams.w * PSOparams.wdamp;
end

out = GlobalBest; % include position, cost, senQuality, mainCost
out.BestCosts = BestCosts; % cost iteration curve
end

