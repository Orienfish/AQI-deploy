%% main.m
clc;
clear;
close all;
warning('off','all')
addpath('./libs/');
addpath('./mlibs/');
addpath('./lldistkm/');
addpath('./gp/');
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

%% get the estimated mean and cov at V
[pm2_5_mean_vd, pm2_5_cov_vd] = gp_predict_knownD(V, D, mean_pm2_5, cov_mat_pm2_5, K_pm2_5);
%bubbleplot_wsize(D(:, 1), D(:, 2), mean_pm2_5, var_pm2_5, 'pm2.5 of D');
%bubbleplot_wsize(V(:, 1), V(:, 2), pm2_5_mean_vd, diag(pm2_5_cov_vd), 'pm2.5 V given D');

%temp_mean_vd = temp_mean_vd / 4 + 180; % weird fix
%bubbleplot_wsize(D(:, 1), D(:, 2), mean_temp, var_temp, 'temp of D');
%bubbleplot_wsize(V(:, 1), V(:, 2), temp_mean_vd, diag(temp_cov_vd), 'temp V given D');

%% setting parameters for algorithms
% setting Quality parameters
Qparams.Xv = V;                         % list of reference locations to 
                                        % predict, [lat lon]
Qparams.cov_vd = pm2_5_cov_vd;          % cov matrix at Xv given pre-deployment D
Qparams.Xd = D;                         % list of predeployment locations
Qparams.mean_d = mean_pm2_5;            % mean value at D
Qparams.cov_d = cov_mat_pm2_5;          % cov matrix at D
Qparams.mean_temp_d = mean_temp;        % mean temperature at D
Qparams.cov_temp_d = cov_mat_temp;      % cov matrix of temperature at D

% set parameters
params.m_A = m_A;                       % number of sensors to deploy
params.Cm = Cm;                         % maintenance cost budget
params.K = K_pm2_5;                     % the fitted RBF kernel function
params.K_temp = K_temp;                 % the fitted RBF kernel function 
                                        % for temperature
params.c = c;                           % position of the sink in [lat lon]
params.R = R;                           % communication range of the sensors in km
params.bound = bound;                   % bound for the area
params.logging = false;                 % logging flag
% parameters of the cost function
params.weights = [0.5 0.4 0.1];         % weights for sensing quality,
                                        % maintenance cost and penalty
params.penalty = 100;                   % penalty for non-connected nodes

%% call the greedy heuristic IDSQ
fprintf('Calling IDSQ...\n');
IDSQparams.alpha = 0.6;             % the weight factor in IDSQ
%resIDSQ = IDSQ(Qparams, params, IDSQparams);
%plot_IDSQ(resIDSQ.Xa, resIDSQ.commMST, c);
%fprintf('IDSQ: senQ: %f mainCost: %f\n', resIDSQ.F, resIDSQ.M);

%% call PSO
fprintf('Calling PSO...\n');
% problem definition
PSOparams.nVar = m_A;                   % number of unknown decision variables
PSOparams.VarSize = [m_A 2]; % matrix size of decision variables
% parameters of PSO
PSOparams.maxIter = 100;                % maximum number of iterations
PSOparams.nPop = 50;                    % populaton size
PSOparams.chi = 0.729;                  % constriction factor
PSOparams.w = PSOparams.chi;            % inertia coefficient
PSOparams.wdamp = 1;                    % damping ratio of inertia coefficient
PSOparams.c1 = 2 * PSOparams.chi;       % personal acceleration coefficient
PSOparams.c2 = 2 * PSOparams.chi;       % social acceleration coefficient

params.Cm = params.Cm - 0.5;            % need some margin

%resPSO = PSO(Qparams, params, PSOparams);

% plot the BestCosts curve
%figure();
%plot(resPSO.BestCosts, 'LineWidth', 2);
%xlabel('Iteration');
%ylabel('Best Cost');

% plot the solution
%[PSOTree, PSOpred] = MST(resPSO.Position, c, R);
%nodesPSO = vertcat(resPSO.Position, c);
%plot_solution(nodesPSO, PSOpred);


%% call ABC
fprintf('Calling ABC...\n');
% problem definition
ABCparams.nVar = m_A;                   % number of unknown decision variables
ABCparams.VarSize = [m_A 2]; % matrix size of decision variables
% parameters of ABC
ABCparams.maxIter = 100;                % maximum number of iterations
ABCparams.nPop = 50;                    % populaton size
ABCparams.nOnlooker = ABCparams.nPop;   % number of onlooker bees
ABCparams.L = round(0.6 * ABCparams.nVar * ABCparams.nPop); 
                                        % Abandonment Limit Parameter (Trial Limit)
ABCparams.a = 1;                        % Acceleration Coefficient Upper Bound
                                        
resABC = ABC(Qparams, params, ABCparams);

% plot the BestCosts curve
figure();
plot(resABC.BestCosts, 'LineWidth', 2);
xlabel('Iteration');
ylabel('Best Cost');

% plot the solution
[ABCTree, ABCpred] = MST(resABC.Position, c, R);
nodesABC = vertcat(resABCA.Position, c);
plot_solution(nodesABC, ABCpred);

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

function plot_solution(nodes, pred)
    figure();
    geoaxes('NextPlot','add');
    % plot connections from commMST
    for i = 1:length(pred)
        if ~isnan(pred(i)) && pred(i) > 0
            geoplot([nodes(i, 1) nodes(pred(i), 1)], ...
                [nodes(i, 2) nodes(pred(i), 2)], 'b-*', 'LineWidth', 2);
        end
    end   
end
