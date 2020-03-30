%% main.m
clc;
clear;
close all;
warning('off','all')
addpath('./libs/');
addpath('./mlibs/');
addpath('./lldistkm/');
diary 'log.txt';

%% settings
% D - pre-deployment
% V - reference locations
% A - deployment plan
m_A = 20;
Cm = 8;
sdate = '2019-01-01 00:00:00 UTC'; % start date of the dataset
edate = '2020-02-20 23:50:00 UTC'; % end date of the dataset
thres = 1e3;                       % a threshold used to filter out outliers
interval = 60 * 10;                % 10 mins = 600 secs
R = 10;                            % communication range of sensors in km

%% pre-process
fprintf('start pre-processing...\n');
% get the mean, var and count of each type of data
f_list = dir('./data/*Primary*.csv'); % use all primary data
dataT_sav = './data/dataT.csv';
tic
if exist(dataT_sav, 'file')
    dataT = readtable(dataT_sav);
else
    % if no file exists yet, do pre-process and save it to file
    dataT = preprocess(f_list, dataT_sav, thres);
end
toc

% get pre-deployment D
D_lat = vertcat(dataT.lat(:));
D_lon = vertcat(dataT.lon(:));
D = horzcat(D_lat, D_lon);

% get the region range, set grid unit and obtain V
bound.latUpper = max(D_lat);
bound.latLower = min(D_lat);
bound.lonUpper = max(D_lon);
bound.lonLower = min(D_lon);
gridUnit = 0.05; % adjustable
c = [(bound.latUpper + bound.latLower) / 2, (bound.lonUpper + bound.lonLower) / 2]; % sink position
n_latV = floor((bound.latUpper - bound.latLower) / gridUnit);
n_lonV = floor((bound.lonUpper - bound.lonLower) / gridUnit);
V_lat = linspace(bound.latLower, bound.latUpper, n_latV);
V_lon = linspace(bound.lonLower, bound.lonUpper, n_lonV);

% fill the grid locations into V
n_V = n_latV * n_lonV; % total number of locations in V
V = zeros(n_V, 2);
for i = 1:n_latV
    for j = 1:n_lonV
        V((i - 1) * n_lonV + j, :) = [V_lat(i) V_lon(j)];
    end
end
fprintf('Generate V with size %d x %d\n', num2str(n_latV), num2str(n_lonV));
%bubbleplot_wsize(V(:, 1), V(:, 2), 1:n_V, 1:n_V, 'order');

% obtain mean vector and covariance matrix and correlation matrix of certain data types
% pm2_5
mean_pm2_5 = vertcat(dataT.pm2_5_avg(:)); % mean
var_pm2_5 = vertcat(dataT.pm2_5_var(:)); % var
cov_mat_pm2_5_sav = './data/cov_mat_pm2_5.csv';
corr_mat_pm2_5_sav = './data/corr_mat_pm2_5.csv';
tic
if exist(cov_mat_pm2_5_sav, 'file') && exist(corr_mat_pm2_5_sav, 'file')
    cov_mat_pm2_5 = readmatrix(cov_mat_pm2_5_sav);
    corr_mat_pm2_5 = readmatrix(corr_mat_pm2_5_sav);
else
[cov_mat_pm2_5, corr_mat_pm2_5] = cov_matrix(f_list, 'pm2_5', ...
    sdate, edate, interval, cov_mat_pm2_5_sav, corr_mat_pm2_5_sav, thres);
end
toc

% temp
mean_temp = vertcat(dataT.temp_avg(:)); % mean
var_temp = vertcat(dataT.temp_var(:)); % var
cov_mat_temp_sav = './data/cov_mat_temp.csv';
corr_mat_temp_sav = './data/corr_mat_temp.csv';
tic
if exist(cov_mat_temp_sav, 'file') && exist(corr_mat_temp_sav, 'file')
    cov_mat_temp = readmatrix(cov_mat_temp_sav);
    corr_mat_temp = readmatrix(corr_mat_temp_sav);
else
[cov_mat_temp, corr_mat_temp] = cov_matrix(f_list, 'temp', ...
    sdate, edate, interval, cov_mat_temp_sav, corr_mat_temp_sav, thres);
end
toc
%figure();
%h = heatmap(cov_mat_temp);

%% fit the RBF kernel for certain data types
fprintf('Fitting the RBF kernel...\n');
K_pm2_5 = fit_kernel(dataT.lat, dataT.lon, cov_mat_pm2_5, 'pm2.5');
K_temp = fit_kernel(dataT.lat, dataT.lon, cov_mat_temp, 'temp');

%% Get the estimated mean and cov at V
[pm2_5_mean_vd, pm2_5_cov_vd] = gp_predict_knownD(V, D, mean_pm2_5, cov_mat_pm2_5, K_pm2_5);
%bubbleplot_wsize(D(:, 1), D(:, 2), mean_pm2_5, var_pm2_5, 'pm2.5 of D');
%bubbleplot_wsize(V(:, 1), V(:, 2), pm2_5_mean_vd, diag(pm2_5_cov_vd), 'pm2.5 V given D');

[temp_mean_vd, temp_cov_vd] = gp_predict_knownD(V, D, mean_temp, cov_mat_temp, K_temp);
%Xv = V;
%Xd = D;
%K = K_temp;
%mean_d = mean_temp;
%cov_d = cov_mat_temp;
%cov_d_inv = inv(cov_d);
% calculate Sigma_VD, Sigma_DV and Sigma_VV
%Sigma_VD = gen_Sigma(Xv, Xd, K);
%Sigma_DV = Sigma_VD';
%Sigma_VV = gen_Sigma(Xv, Xv, K);

% calculate mean vector and covariance matrix
%temp_mean_vd = Sigma_VD * cov_d_inv * mean_d;
%temp_cov_vd = Sigma_VV - Sigma_VD * cov_d_inv * Sigma_DV;
temp_mean_vd = temp_mean_vd / 4 + 180; % weird fix
%bubbleplot_wsize(D(:, 1), D(:, 2), mean_temp, var_temp, 'temp of D');
%bubbleplot_wsize(V(:, 1), V(:, 2), temp_mean_vd, diag(temp_cov_vd), 'temp V given D');

% Generate a random A
%A = zeros(m_A, 2);
%pd_lat = makedist('Uniform', 'lower', bound.latLower, 'upper', bound.latUpper);
%pd_lon = makedist('Uniform', 'lower', bound.lonLower, 'upper', bound.lonUpper);
%for i = 1:m_A
%    A(i, :) = [random(pd_lat) random(pd_lon)];
%end
% plot V and A on one map
%sizedata = ones(m_A + n_V, 1); % uniform size
%colordata = categorical(vertcat(zeros(m_A, 1), ones(n_V, 1))); % two categories
%bubbleplot_wcolor(vertcat(A(:, 1), V(:, 1)), vertcat(A(:, 2), V(:, 2)), ...
%    sizedata, colordata, 'V and A');
% get the estimated mean and cov at A
%[mean_ad, cov_ad] = gp_predict_knownA(A, D, mean_pm2_5, cov_mat_pm2_5, K);
%bubbleplot_wsize(A(:, 1), A(:, 2), mean_ad, diag(cov_ad), 'A given D');

% Evaluate uncertainty or conditional entropy
%condEntropy = cond_entropy(V, A, K);
%condEntropy_d = cond_entropy_d(V, cov_vd, A, cov_ad, K);
%fprintf('conditional entropy w/o predeployment: %f\n', condEntropy);
%fprintf('conditional entropy w/ predeployment: %f\n', condEntropy_d);
%senseQuality = sense_quality(V, cov_vd, A, cov_ad, K);
%fprintf('sensing quality w/ predeployment: %f\n', senseQuality);

%% Call the greedy heuristic IDSQ
fprintf('Calling IDSQ...\n');
Tv_cel = fah2cel(temp_mean_vd);
params.m_A = m_A;
params.Cm = Cm;
params.K = K_pm2_5;
params.c = c;
params.R = R;
params.IDSQ_alpha = 0.6;
params.logging = false;
%resIDSQ = IDSQ(V, pm2_5_cov_vd, Tv_cel, params);
%plot_IDSQ(resIDSQ.Xa, resIDSQ.commMST, c);

%% call PSO
fprintf('Calling PSO...\n');
% problem definition
nVar = m_A;             % number of unknown decision variables
VarSize = [nVar 2];     % matrix size of decision variables

% parameters of PSO
maxIter = 100;          % maximum number of iterations
nPop = 50;              % populaton size
w = 1;                  % inertia coefficient
wdamp = 0.99;           % damping ratio of inertia coefficient
c1 = 2;                 % personal acceleration coefficient
c2 = 2;                 % social acceleration coefficient
% parameters of the cost function
params.weights = [10 1 0.2];
params.penalty = 100;

% initialization
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
particle = repmat(empty_particle, nPop, 1);

% initialize global best
GlobalBest.Cost = inf;
GlobalBest.senQuality = inf;
GlobalBest.mainCost = inf;

% initialize population members
for i=1:nPop
    % generate random solution
    particle(i).Position(:, 1) = unifrnd(bound.latLower, bound.latUpper, nVar, 1);
    particle(i).Position(:, 2) = unifrnd(bound.lonLower, bound.lonUpper, nVar, 1);
    
    % initialize velocity
    particle(i).Velocity = zeros(VarSize);
    
    % cost evaluation
    [pm2_5_mean_ad, pm2_5_cov_ad] = gp_predict_knownD( ...
        particle(i).Position, D, mean_pm2_5, cov_mat_pm2_5, K_pm2_5);
    [temp_mean_ad, temp_cov_ad] = gp_predict_knownD( ...
        particle(i).Position, D, mean_temp, cov_mat_temp, K_temp);
    temp_mean_ad = temp_mean_ad / 4 + 180; % weird fix
    
    % setting Quality parameters
    Qparams.Xv = V;
    Qparams.cov_vd = pm2_5_cov_vd;
    Qparams.Xa = particle(i).Position;
    Qparams.Ta = fah2cel(temp_mean_ad);
    Qparams.cov_ad = pm2_5_cov_ad;
    [res, commMST, predMST] = costFunction(Qparams, params);
    
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
BestCosts = zeros(maxIter, 1);

% main loop of PSO
for it = 1:maxIter
    for i = 1:nPop
        % update velocity
        particle(i).Velocity = particle(i).Velocity ...
            + c1*rand(VarSize).*(particle(i).Best.Position - particle(i).Position) ...
            + c2*rand(VarSize).*(GlobalBest.Position - particle(i).Position);
        
        % update position
        particle(i).position = particle(i).Position + particle(i).Velocity;
        
        % cost evaluation
        [pm2_5_mean_ad, pm2_5_cov_ad] = gp_predict_knownD( ...
            particle(i).Position, D, mean_pm2_5, cov_mat_pm2_5, K_pm2_5);
        [temp_mean_ad, temp_cov_ad] = gp_predict_knownD( ...
            particle(i).Position, D, mean_temp, cov_mat_temp, K_temp);
        temp_mean_ad = temp_mean_ad / 4 + 180; % weird fix
        
        % Qparams.Xv = V;                 % use the same value as init
        % Qparams.cov_vd = pm2_5_cov_vd;  % use the same value as init
        Qparams.Xa = particle(i).Position;
        Qparams.Ta = fah2cel(temp_mean_ad);
        Qparams.cov_ad = pm2_5_cov_ad;
        [cost, commMST, predMST] = costFunction(Qparams, params);

        particle(i).Cost = cost;
        particle(i).senQuality = res.F;
        particle(i).mainCost = res.M;
        
        % update personal best
        if particle(i).Cost < particle(i).Best.Cost
            particle(i).Best.Position = paticle(i).Position;
            particle(i).Best.Cost = particle(i).Cost;
            particle(i).Best.senQuality = particle(i).senQuality;
            particle(i).Best.mainCost = particle(i).mainCost;
            
            % update the global best
            if particle(i).Best.Cost < GlobalBest.Cost
                GlobalBest = particle(i).Best;
            end
        end
    end
    
    % store the best cost value
    BestCosts(it) = GlobalBest.Cost;
    
    % display iteration information
    fprintf('Iteration %d: Best Cost: %f senQ: %f mainCost: %f\n', ...
        it, GlobalBest.Cost, GlobalBest.senQuality, GlobalBest.mainCost);
    disp([])
    
    % damping inertia coefficient
    w = w * wdamp;
end

figure;
plot(BestCosts, 'LineWidth', 2);
xlabel('Iteration');
ylabel('Best Cost');
% plot_solution(A, predMST, c);


%% plot functions
function bubbleplot(lat, lon, title)
    % plot the locations
    figure('Position', [0 0 1000 800]);
    geobubble(lat, lon, 'Title', title);
    %geobasemap streets-light; % set base map style
end

function bubbleplot_wsize(lat, lon, mean, var, title)
    % plot the locations
    figure('Position', [0 0 2000 800]);
    subplot(1, 2, 1);
    geobubble(lat, lon, mean, 'Title', title);
    subplot(1, 2, 2);
    geobubble(lat, lon, var, 'Title', title);
    %geobasemap streets-light; % set base map style
end

function bubbleplot_wcolor(lat, lon, sizedata, colordata, title)
    % plot the locations
    figure('Position', [0 0 1000 800]);
    geobubble(lat, lon, sizedata, colordata, 'Title', title);
    %geobasemap streets-light; % set base map style
end

function plot_IDSQ(A, commMST, c)
    figure();
    geoaxes('NextPlot','add');
    % generate list of nodes including the sink
    nodes = vertcat(A, c);
    % plot connections from commMST
    for i = 1:size(commMST, 1)
        for j = 1:size(commMST, 2)
            if ~isnan(commMST(i, j))
                geoplot([nodes(i, 1) nodes(j, 1)], ...
                    [nodes(i, 2) nodes(j, 2)], 'b-*', 'LineWidth', 2);
            end
        end
    end   
end

function plot_solution(A, predMST, c)
    figure();
    geoaxes('NextPlot','add');
    % generate list of nodes including the sink
    nodes = vertcat(A, c);
    % plot connections from commMST
    for i = 1:length(predMST)
        if ~isnan(predMST(i)) && predMST(i) > 0
            geoplot([nodes(i, 1) nodes(predMST(i), 1)], ...
                [nodes(i, 2) nodes(predMST(i), 2)], 'b-*', 'LineWidth', 2);
        end
    end   
end
