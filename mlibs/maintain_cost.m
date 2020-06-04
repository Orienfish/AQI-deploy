function [out] = maintain_cost(Xa, Ta, conn_idx, commMST, predMST, logging)
%% Compute the maintenance cost from a given deployment.
%  Note deployment Xa is a complete list. Not all nodes are connected.
%
% Args:
%   Xa: a given deployment list, [lat lon]
%   Ta: average temperature estimation at Xa in Celsius 
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

% settings for battery
cap_bat = 20000;   % initial battery capacity in mAh
dt_bat_h = 1;      % time resolution of battery in hours
c_bat = 10;        % cost to replace battery

% setting for circuit
c_node = 100;     % cost to replace node
C = 0;            % total maintenance cost

% init return lifetime list
out.batlife = zeros(size(Xa, 1), 1);
out.cirlife = zeros(size(Xa, 1), 1);

% get children cnt of each node from MST
child_cnt = get_child_cnt(predMST);

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
        
        % estimate power in W
        [stbPwr, stbTc] = stbPower(params, Ta(i));

        % estimate battery lifetime in days
        I_mA = stbPwr * 1000 / params.Vdd; % convert from W to mW then calculate average current draw
        batlife_ratio = bat_ratio(cap_bat, Ta(i), I_mA, dt_bat_h);

        % estimate circuit lifetime in days
        cirlife_ratio = mttf_ratio(stbTc);
        
        % update total maintenance cost
        nodeC = c_bat / batlife_ratio + c_node / cirlife_ratio;
        C = C + nodeC;
        
        % update output list
        out.batlife(i) = batlife_day;
        out.cirlife(i) = cirlife_day;

        if logging
            fprintf('  node %d amb temp: %f avg pwr: %f core temp: %f\n', ...
                i, Ta(i), stbPwr, stbTc);
            fprintf('  bat life: %f cir life: %f main cost: %f\n', ...
                batlife_day, cirlife_day, nodeC);
        end
    end
end

out.C = C;
end

