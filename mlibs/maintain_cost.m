function [out] = maintain_cost(Xa, Ta, Ta_v, conn_idx, commMST, predMST, logging)
%% Compute the maintenance cost from a given deployment.
%  Note deployment Xa is a complete list. Not all nodes are connected.
%
% Args:
%   Xa: a given deployment list, [lat lon]
%   Ta: average temperature estimation at Xa in Celsius
%   Ta_v: variance of ambient temperature at Xa
%   conn_idx: list of logical variables showing connection
%   commMST: generated feasible MST for connnection
%   predMST: the predecessor of the node
%   logging: logging flag
%
% Return:
%   out.C: the maintenance cost
%   out.batlife: list of battery lifetime in days
%   out.cirlife: list of circuit lifetime in days

% settings for sensor workloads
params.Pto = 0.52;       % 520mW transmission power baseline
params.Btx = 125;        % 1kbps = 125B/s bandwidth
params.Brx = 125;
params.Ltx_init = 100;   % 100B packet length
params.Lrx_init = 100;           
params.Prx = 0.2;        % 200mW receiving power
params.Psen = 0.2;       % 200mW senisng power
params.tsen = 0.3;       % 300ms sensing time
params.T = 10;           % 10s sampling frequency
params.f = 300e6;        % 300MHz clock frequency
params.Vdd = 3.3;        % 3.3v supply voltage
params.Iref = 50;        % 50mA reference current draw
params.Pref = params.Vdd * params.Iref / 1000; % reference power in W
params.nbins = 10;       % number of bins to deal with temperature variation

% settings for battery
params.cap_bat = 2000;   % initial battery capacity in mAh
params.dt_bat_h = 1;     % time resolution of battery in hours
params.c_bat = 10;       % cost to replace battery

% setting for circuit
params.c_node = 100;     % cost to replace node

% init return list
out.C = 0.0;             % total maintenance cost
out.batlife = zeros(size(Xa, 1), 1);
out.cirlife = zeros(size(Xa, 1), 1);

% get children cnt of each node from MST
child_cnt = get_child_cnt(predMST);

% temporary fix
Ta_v = Ta_v * 5;

% iterative through every node in Xa
for i = 1:size(Xa, 1)
    if conn_idx(i) % only process those nodes that are connected
        txDist_km = commMST(predMST(i), i); % get the distance to predecessor
        % txDist_m = txDist_km / 1000;  % convert to meters
        if isinf(txDist_km)
            error('Inf commMST!'); % security check
        end
        params.dtx = txDist_km;
        params.Ltx = params.Ltx_init * (child_cnt(i) + 1);
        params.Lrx = params.Lrx_init * child_cnt(i);
        % fprintf('child cnt: %d Ltx: %d Lrx: %d\n', child_cnt, params.Ltx, params.Lrx);
        
        % deal with temperature deviation
        % simulate temperature distribution and count the prob. in each bin
        curCenters = linspace(Ta(i)-3*Ta_v(i), Ta(i)+3*Ta_v(i), params.nbins);
        T_cdf = [normcdf(curCenters, Ta(i), Ta_v(i)), 1.0];
        curProb = T_cdf(2:end) - T_cdf(1:end-1);
        
        [nodeC, batlife_ratio, cirlife_ratio] = maintain_node(params, ...
            curCenters, curProb);
        % update output list
        out.C = out.C + nodeC;
        out.batlife(i) = batlife_ratio;
        out.cirlife(i) = cirlife_ratio;

        if logging
            fprintf('  node %d amb temp: %f avg pwr: %f core temp: %f\n', ...
                i, Ta(i), stbPwr, stbTc);
            fprintf('  bat life: %f cir life: %f main cost: %f\n', ...
                batlife_ratio, cirlife_ratio, nodeC);
        end
        
        % for comparison, the following line compute the original 
        % maintenance cost only using average temperature
        % estimate power in W
        %[stbPwr, stbTc] = stbPower(params, Ta(i));

        % estimate battery lifetime in ratio
        % convert from W to mW then calculate average current draw
        %I_mA = stbPwr * 1000 / params.Vdd;
        %batlife_cur = bat_ratio(cap_bat, Ta(i), I_mA, dt_bat_h);

        % estimate circuit lifetime in ratio
        %cirlife_cur = mttf_ratio(stbTc);

        % update total maintenance cost
        %nodeCur = c_bat / batlife_cur + c_node / cirlife_cur;
        %fprintf('    original maintenance cost: %f\n', nodeCur);
    end
end
end

