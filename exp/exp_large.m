function exp_large(target, run, Q)
%% Run simulation on the large dataset for the given target
%
% Args:
%   run: experimental parameters setting
%   Q: pre-determined sensing quality quota
warning('off','all')
addpath('../libs/');
addpath('../mlibs/');
addpath('../lldistkm/');
addpath('../gp/');
addpath('../');

%% settings
% D - pre-deployment
% V - reference locations
% A - deployment plan
m_A = 16;                          % number of sensors to place
%Q = 10.0;                        % sensing quality quota
R = 10;                            % communication range of sensors in km
sdate = '2019-01-01 00:00:00 UTC'; % start date of the dataset
edate = '2020-04-01 23:00:00 UTC'; % end date of the dataset
thres = 1e3;                       % a threshold used to filter out outliers
interval = 60 * 60;                % 60 mins = 3600 secs

% boolean variables deciding whether to run each algorithm
%run.IDSQ = false;
%run.pSPIEL = true;
%run.PSO = true;
%run.ABC = true;
%run.iter = 10;

%% pre-process
fprintf('start pre-processing...\n');
% get the mean, var and count of each type of data
f_list = dir('../dataL/*Primary*.csv'); % use all primary data
dataT_sav = '../dataL/dataT.csv';
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
fprintf('Generate V with size %d x %d\n', n_latV, n_lonV);
%bubbleplot_wsize(V(:, 1), V(:, 2), 1:n_V, 'order');

%% obtain mean vector and covariance matrix and correlation matrix of certain data types
% We need temperature data anyway
mean_temp = vertcat(dataT.temp_avg(:)); % mean
var_temp = vertcat(dataT.temp_var(:)); % var
cov_mat_temp_sav = '../dataL/cov_mat_temp.csv';
corr_mat_temp_sav = '../dataL/corr_mat_temp.csv';
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

% pm2_5
if strcmp(target, 'pm2_5')
mean_target = vertcat(dataT.pm2_5_avg(:)); % mean
var_target = vertcat(dataT.pm2_5_var(:)); % var
cov_mat_target_sav = '../dataL/cov_mat_pm2_5.csv';
corr_mat_target_sav = '../dataL/corr_mat_pm2_5.csv';
tic
if exist(cov_mat_target_sav, 'file') && exist(corr_mat_target_sav, 'file')
    cov_mat_target = readmatrix(cov_mat_target_sav);
    corr_mat_target = readmatrix(corr_mat_target_sav);
else
[cov_mat_target, corr_mat_target] = cov_matrix(f_list, 'pm2_5', ...
    sdate, edate, interval, cov_mat_target_sav, corr_mat_target_sav, thres);
end
toc
end

% temp
if strcmp(target, 'temp')
mean_target = mean_temp; % mean
var_target = var_temp; % var
cov_mat_target = cov_mat_temp;
corr_mat_target = corr_mat_temp;
end

% pm1
if strcmp(target, 'pm1')
mean_target = vertcat(dataT.pm1_avg(:)); % mean
var_target = vertcat(dataT.pm1_var(:)); % var
cov_mat_target_sav = '../dataL/cov_mat_pm1.csv';
corr_mat_target_sav = '../dataL/corr_mat_pm1.csv';
tic
if exist(cov_mat_target_sav, 'file') && exist(corr_mat_target_sav, 'file')
    cov_mat_target = readmatrix(cov_mat_target_sav);
    corr_mat_target = readmatrix(corr_mat_target_sav);
else
[cov_mat_target, corr_mat_target] = cov_matrix(f_list, 'pm1', ...
    sdate, edate, interval, cov_mat_target_sav, corr_mat_target_sav, thres);
end
toc
end

% pm10
if strcmp(target, 'pm10')
mean_target = vertcat(dataT.pm10_avg(:)); % mean
var_target = vertcat(dataT.pm10_var(:)); % var
cov_mat_target_sav = '../dataL/cov_mat_pm10.csv';
corr_mat_target_sav = '../dataL/corr_mat_pm10.csv';
tic
if exist(cov_mat_target_sav, 'file') && exist(corr_mat_target_sav, 'file')
    cov_mat_target = readmatrix(cov_mat_target_sav);
    corr_mat_target = readmatrix(corr_mat_target_sav);
else
[cov_mat_target, corr_mat_target] = cov_matrix(f_list, 'pm10', ...
    sdate, edate, interval, cov_mat_target_sav, corr_mat_target_sav, thres);
end
toc
end

% humid
if strcmp(target, 'humid')
mean_target = vertcat(dataT.humid_avg(:)); % mean
var_target = vertcat(dataT.humid_var(:)); % var
cov_mat_target_sav = '../dataL/cov_mat_humid.csv';
corr_mat_target_sav = '../dataL/corr_mat_humid.csv';
tic
if exist(cov_mat_target_sav, 'file') && exist(corr_mat_target_sav, 'file')
    cov_mat_target = readmatrix(cov_mat_target_sav);
    corr_mat_target = readmatrix(corr_mat_target_sav);
else
[cov_mat_target, corr_mat_target] = cov_matrix(f_list, 'humid', ...
    sdate, edate, interval, cov_mat_target_sav, corr_mat_target_sav, thres);
end
toc
end


%% fit the RBF kernel for certain data types
fprintf('Fitting the RBF kernel...\n');
K_target = fit_kernel(dataT.lat, dataT.lon, cov_mat_target, target);
K_temp = fit_kernel(dataT.lat, dataT.lon, cov_mat_temp, 'temp');
% shifting temp matrix to standardized one at predeployment locations
cov_mat_temp = gen_Sigma(D, D, K_temp);
fix matrix

%% get the estimated mean and cov at V for plotting
if run.debugPlot
    [target_mean_vd, target_cov_vd] = gp_predict_knownD(V, D, mean_target, ...
        cov_mat_target, K_target);
    bubbleplot_wsize(D(:, 1), D(:, 2), mean_target, 'mean of target at D');
    bubbleplot_wsize(D(:, 1), D(:, 2), var_target, 'variance of target at D');
    bubbleplot_wsize(V(:, 1), V(:, 2), target_mean_vd, 'mean of target at V given D');
    bubbleplot_wsize(V(:, 1), V(:, 2), diag(target_cov_vd), ...
        'variance of target at V given D');

    [temp_mean_vd, temp_cov_vd] = gp_predict_knownD(V, D, mean_temp, ...
        cov_mat_temp, K_temp);
    temp_mean_vd = fah2cel(temp_mean_vd); % convert to Celsius
    % for plotting distribution
    bubbleplot_wsize(D(:, 1), D(:, 2), mean_temp, 'mean of temp at D');
    bubbleplot_wsize(D(:, 1), D(:, 2), var_temp, 'variance of temp at D');
    bubbleplot_wsize(V(:, 1), V(:, 2), temp_mean_vd, 'mean of temp at V given D');
    bubbleplot_wsize(V(:, 1), V(:, 2), diag(temp_cov_vd), ...
        'variance of temp at V given D');
    % plot heatmap of temperature
    temp_mean_vd = flipud(reshape(temp_mean_vd, [n_lonV, n_latV])');
    figure;
    h = heatmap(round(V_lon*100)/100, round(V_lat*100)/100, ...
        temp_mean_vd, 'Colormap', flipud(autumn), 'CellLabelColor','none', ...
        'XLabel','Longitude', 'YLabel', 'Latitude', 'FontSize', 16);
end

%% setting parameters for algorithms
% setting Quality parameters
Qparams.Xv = V;                         % list of reference locations to 
                                        % predict, [lat lon]
Qparams.cov_vd = gen_Sigma(V, V, K_target); % cov matrix at Xv given pre-deployment D
Qparams.Xd = D;                         % list of predeployment locations
Qparams.mean_d = mean_target;           % mean value at D
Qparams.cov_d = cov_mat_target;         % cov matrix at D
Qparams.mean_temp_d = mean_temp;        % mean temperature at D
Qparams.cov_temp_d = cov_mat_temp;      % cov matrix of temperature at D

% set parameters
params.n_V = n_V;                       % number of reference locations
params.m_A = m_A;                       % number of sensors to deploy
params.Q = Q;                           % sensing quality quota
params.K = K_target;                    % the fitted RBF kernel function
params.K_temp = K_temp;                 % the fitted RBF kernel function 
                                        % for temperature
params.c = c;                           % position of the sink in [lat lon]
params.R = R;                           % communication range of the sensors in km
params.bound = bound;                   % bound for the area
params.logging = false;                 % logging flag
% parameters of the cost function
params.weights = [0.5 0.4 0.1];         % weights for maintenance cost,
                                        % sensing quality and penalty
params.penalty = 100;                   % penalty for non-connected nodes

%% call the greedy heuristic IDSQ
if run.IDSQ
    for it = 1:run.iter
        fprintf('Calling IDSQ...\n');
        IDSQparams.alpha = 0.6;             % the weight factor in IDSQ
        resIDSQ = IDSQ(Qparams, params, IDSQparams);
        plot_IDSQ(resIDSQ.Xa, resIDSQ.commMST, c);
        fprintf('IDSQ: senQ: %f mainCost: %f\n', resIDSQ.F, resIDSQ.M);
    end
end

%% call pSPIEL
if run.pSPIEL
    for it = 1:run.iter
        fprintf('Calling pSPIEL...\n');
        tic
        respSPIEL = pSPIEL(Qparams, params);
        toc
        fprintf('pSPIEL: # of nodes: %d senQ: %f mainCost: %f\n', ...
            sum(respSPIEL.connected), respSPIEL.F, respSPIEL.M.C);
        nodespSPIEL = vertcat(respSPIEL.Position, c);
        %plot_solution(nodespSPIEL, respSPIEL.pred);
        %bubbleplot_wsize(respSPIEL.Position(:, 1), respSPIEL.Position(:, 2), ...
        %    respSPIEL.M.batlife, '');
        %bubbleplot_wsize(respSPIEL.Position(:, 1), respSPIEL.Position(:, 2), ...
        %    respSPIEL.M.cirlife, '');
        
        % logging
        bat_str = '';
        cir_str = '';
        for idx = 1:params.n_V
            bat_str = sprintf('%s%.4f,', bat_str, respSPIEL.M.batlife(idx));
            cir_str = sprintf('%s%.4f,', cir_str, respSPIEL.M.cirlife(idx));
        end
        str = sprintf('%d %f %f', sum(respSPIEL.connected), respSPIEL.F, ...
            respSPIEL.M.C);
        str = sprintf('%s\n%s\n%s\n', str, bat_str, cir_str);
        filename = sprintf('pSPIEL_%s_%d.txt', target, Q);
        log(filename, str);
    end
end

%% call PSO
if run.PSO
    for it = 1:run.iter
    fprintf('Calling PSO...\n');
        % problem definition
        PSOparams.nVar = m_A;                   % number of unknown decision variables
        PSOparams.VarSize = [m_A 2]; % matrix size of decision variables
        % parameters of PSO
        PSOparams.maxIter = 50;                % maximum number of iterations
        PSOparams.nPop = 50;                    % populaton size
        PSOparams.chi = 0.729;                  % constriction factor
        PSOparams.w = PSOparams.chi;            % inertia coefficient
        PSOparams.wdamp = 1;                    % damping ratio of inertia coefficient
        PSOparams.c1 = 2 * PSOparams.chi;       % personal acceleration coefficient
        PSOparams.c2 = 2 * PSOparams.chi;       % social acceleration coefficient

        tic
        resPSO = PSO(Qparams, params, PSOparams);
        toc

        % plot the BestCosts curve
        %figure();
        %plot(resPSO.BestCosts, 'LineWidth', 2);
        %xlabel('Iteration');
        %ylabel('Best Cost');

        % plot the solution
        [PSO_G, PSOpred] = MST(resPSO.Position, c, R);
        nodesPSO = vertcat(resPSO.Position, c);
        %plot_solution(nodesPSO, PSOpred);

        % plot lifetime of each node
        connected = ~isnan(PSOpred(1:params.m_A)); % a logical array of connected sensors
        [temp_mean_ad, temp_cov_ad] = gp_predict_knownD( ...
            resPSO.Position, Qparams.Xd, Qparams.mean_temp_d, ...
            Qparams.cov_temp_d, params.K_temp);
        Qparams.Xa = resPSO.Position;
        Qparams.Ta = fah2cel(temp_mean_ad);

        % calculate the maintenance cost of connected sensors
        M = maintain_cost(Qparams.Xa, Qparams.Ta, connected, PSO_G, PSOpred, ...
            params.logging);
        %bubbleplot_wsize(Qparams.Xa(:, 1), Qparams.Xa(:, 2), M.batlife, '');
        %bubbleplot_wsize(Qparams.Xa(:, 1), Qparams.Xa(:, 2), M.cirlife, '');
        
        % logging
        bat_str = '';
        cir_str = '';
        for idx = 1:params.m_A
            bat_str = sprintf('%s%.4f,', bat_str, M.batlife(idx));
            cir_str = sprintf('%s%.4f,', cir_str, M.cirlife(idx));
        end
        str = sprintf('%d %f %f', sum(connected), resPSO.senQuality, ...
            resPSO.mainCost);
        str = sprintf('%s\n%s\n%s\n', str, bat_str, cir_str);
        filename = sprintf('PSO_%s_%d.txt', target, Q);
        log(filename, str);
    end
end

%% call ABC
if run.ABC
    for it = 1:run.iter
        fprintf('Calling ABC...\n');
        % problem definition
        ABCparams.nVar = m_A;                   % number of unknown decision variables
        ABCparams.VarSize = [m_A 2]; % matrix size of decision variables
        % parameters of ABC
        ABCparams.maxIter = 50;                % maximum number of iterations
        ABCparams.nPop = 50;                    % populaton size
        ABCparams.nOnlooker = ABCparams.nPop;   % number of onlooker bees
        ABCparams.L = round(0.4 * ABCparams.nVar * ABCparams.nPop); 
                                                % Abandonment Limit Parameter (Trial Limit)
        ABCparams.a = 0.4;                      % Acceleration Coefficient Upper Bound

        tic
        resABC = ABC(Qparams, params, ABCparams);
        toc

        % plot the BestCosts curve
        %figure();
        %plot(resABC.BestCosts, 'LineWidth', 2);
        %xlabel('Iteration');
        %ylabel('Best Cost');

        % plot the solution
        [ABC_G, ABCpred] = MST(resABC.Position, c, R);
        nodesABC = vertcat(resABC.Position, c);
        %plot_solution(nodesABC, ABCpred);

        % plot lifetime of each node
        connected = ~isnan(ABCpred(1:params.m_A)); % a logical array of connected sensors
        [temp_mean_ad, temp_cov_ad] = gp_predict_knownD( ...
            resABC.Position, Qparams.Xd, Qparams.mean_temp_d, ...
            Qparams.cov_temp_d, params.K_temp);
        Qparams.Xa = resABC.Position;
        Qparams.Ta = fah2cel(temp_mean_ad);

        % calculate the maintenance cost of connected sensors
        M = maintain_cost(Qparams.Xa, Qparams.Ta, connected, ABC_G, ABCpred, ...
            params.logging);
        %bubbleplot_wsize(Qparams.Xa(:, 1), Qparams.Xa(:, 2), M.batlife, '');
        %bubbleplot_wsize(Qparams.Xa(:, 1), Qparams.Xa(:, 2), M.cirlife, '');
        
        % logging
        bat_str = '';
        cir_str = '';
        for idx = 1:params.m_A
            bat_str = sprintf('%s%.4f,', bat_str, M.batlife(idx));
            cir_str = sprintf('%s%.4f,', cir_str, M.cirlife(idx));
        end
        str = sprintf('%d %f %f', sum(connected), resABC.senQuality, ...
            resABC.mainCost);
        str = sprintf('%s\n%s\n%s\n', str, bat_str, cir_str);
        filename = sprintf('ABC_%s_%d.txt', target, Q);
        log(filename, str);
    end
end

end

%% plot functions
function bubbleplot(lat, lon, title)
    % plot the locations
    figure('Position', [0 0 1000 800]);
    geobubble(lat, lon, 'Title', title);
    ax = gca; % get current axes
    ax.FontSize = 16;
    %geobasemap streets-light; % set base map style
end

function bubbleplot_wsize(lat, lon, sizedata, title)
    % plot the locations
    figure;
    geobubble(lat, lon, sizedata, 'Title', title);
    ax = gca; % get current axes
    ax.FontSize = 16;
    %geobasemap streets-light; % set base map style
end

function plot_IDSQ(A, commMST, c)
    figure;
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
    ax = gca; % get current axes
    ax.FontSize = 16;
end

function plot_solution(nodes, pred)
    figure;
    geoaxes('NextPlot','add');
    % plot connections from commMST
    for i = 1:length(pred)
        if ~isnan(pred(i)) && pred(i) > 0
            geoplot([nodes(i, 1) nodes(pred(i), 1)], ...
                [nodes(i, 2) nodes(pred(i), 2)], 'b-*', 'LineWidth', 2);
        end
    end
    ax = gca; % get current axes
    ax.FontSize = 16;   
end

%% log function
function log(filename, str)
    fileID = fopen(filename, 'a+');
    fprintf(fileID, '%s', str);
    fclose(fileID);
end