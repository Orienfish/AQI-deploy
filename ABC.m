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
%   ABCparams.nVar: number of unknown decision variables
%   ABCparams.VarSize: matrix size of decision variables
%   ABCparams.maxIter: maximum number of iterations
%   ABCparams.nPop: population size
%   ABCparams.chi: constriction factor
%   ABCparams.w: inertia coefficient
%   ABCparams.wdamp: damping ratio of inertia coefficient
%   ABCparams.c1: personal acceleration coefficient
%   ABCparams.c2: social acceleration coefficient
%   ABCparams.bound: PSOparams.bound for the area
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


%% ABC Settings
MaxIt=200;              % Maximum Number of Iterations
nPop=100;               % Population Size (Colony Size)
nOnlooker=nPop;         % Number of Onlooker Bees
L=round(0.6*nVar*nPop); % Abandonment Limit Parameter (Trial Limit)
a=1;                    % Acceleration Coefficient Upper Bound
%% Initialization
% Empty Bee Structure
empty_bee.Position=[];
empty_bee.Cost=[];
% Initialize Population Array
pop=repmat(empty_bee,nPop,1);
% Initialize Best Solution Ever Found
BestSol.Cost=inf;
% Create Initial Population
for i=1:nPop
    pop(i).Position=unifrnd(VarMin,VarMax,VarSize);
    pop(i).Cost=CostFunction(pop(i).Position);
    if pop(i).Cost<=BestSol.Cost
        BestSol=pop(i);
    end
end
% Abandonment Counter
C=zeros(nPop,1);
% Array to Hold Best Cost Values
BestCost=zeros(MaxIt,1);
%% ABC Main Loop
for it=1:MaxIt
    
    % Recruited Bees
    for i=1:nPop
        
        % Choose k randomly, not equal to i
        K=[1:i-1 i+1:nPop];
        k=K(randi([1 numel(K)]));
        
        % Define Acceleration Coeff.
        phi=a*unifrnd(-1,+1,VarSize);
        
        % New Bee Position
        newbee.Position=pop(i).Position+phi.*(pop(i).Position-pop(k).Position);
        
        % Evaluation
        newbee.Cost=CostFunction(newbee.Position);
        
        % Comparision
        if newbee.Cost<=pop(i).Cost
            pop(i)=newbee;
        else
            C(i)=C(i)+1;
        end
        
    end
    
    % Calculate Fitness Values and Selection Probabilities
    F=zeros(nPop,1);
    MeanCost = mean([pop.Cost]);
    for i=1:nPop
        F(i) = exp(-pop(i).Cost/MeanCost); % Convert Cost to Fitness
    end
    P=F/sum(F);
    
    % Onlooker Bees
    for m=1:nOnlooker
        
        % Select Source Site
        i=RouletteWheelSelection(P);
        
        % Choose k randomly, not equal to i
        K=[1:i-1 i+1:nPop];
        k=K(randi([1 numel(K)]));
        
        % Define Acceleration Coeff.
        phi=a*unifrnd(-1,+1,VarSize);
        
        % New Bee Position
        newbee.Position=pop(i).Position+phi.*(pop(i).Position-pop(k).Position);
        
        % Evaluation
        newbee.Cost=CostFunction(newbee.Position);
        
        % Comparision
        if newbee.Cost<=pop(i).Cost
            pop(i)=newbee;
        else
            C(i)=C(i)+1;
        end
        
    end
    
    % Scout Bees
    for i=1:nPop
        if C(i)>=L
            pop(i).Position=unifrnd(VarMin,VarMax,VarSize);
            pop(i).Cost=CostFunction(pop(i).Position);
            C(i)=0;
        end
    end
    
    % Update Best Solution Ever Found
    for i=1:nPop
        if pop(i).Cost<=BestSol.Cost
            BestSol=pop(i);
        end
    end
    
    % Store Best Cost Ever Found
    BestCost(it)=BestSol.Cost;
    
    % Display Iteration Information
    disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCost(it))]);
    
end
end

