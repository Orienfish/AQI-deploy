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
params.Btx = 1250;       % 10kbps = 1250B/s bandwidth
params.Brx = 1250;
params.Ltx_init = 1e3;   % 1kB packet length
params.Lrx_init = 1e3;           
params.Prx = 0.2;        % 200mW receiving power
params.Psen = 0.2;       % 200mW senisng power
params.tsen = 0.3;       % 300ms sensing time
params.T = 10;           % 10s sampling frequency
params.f = 300e6;        % 300MHz clock frequency
params.Vdd = 3.3;        % 3.3v supply voltage

% settings for battery
cap_bat = 20000;   % initial battery capacity in mAh
dt_bat_h = 1;     % time resolution of battery in hours
c_bat = 10;        % cost to replace battery

% setting for circuit
c_node = 100;     % cost to replace node
C = 0;            % total maintenance cost

% init return lifetime list
out.batlife = zeros(size(Xa, 1), 1);
out.cirlife = zeros(size(Xa, 1), 1);

% iterative through every node in Xa
for i = 1:size(Xa, 1)
    if conn_idx(i) % only process those nodes that are connected
        txDist_km = commMST(predMST(i), i); % get the distance to predecessor
        % txDist_m = txDist_km / 1000;  % convert to meters
        if isinf(txDist_km)
            error('Inf commMST!'); % security check
        end
        child_cnt = sum(predMST == i); % get child cnt of node i
        params.dtx = txDist_km;
        params.Ltx = params.Ltx_init * (child_cnt + 1);
        params.Lrx = params.Lrx_init * child_cnt;
        % fprintf('child cnt: %d Ltx: %d Lrx: %d\n', child_cnt, params.Ltx, params.Lrx);
        
        % estimate power in W
        [stbPwr, stbTc] = stbPower(params, Ta(i));

        % estimate battery lifetime in days
        I_mA = stbPwr * 1000 / params.Vdd; % convert from W to mW then calculate average current draw
        batlife_h = bat_lifetime(cap_bat, Ta(i), I_mA, dt_bat_h);
        batlife_day = batlife_h / 24;

        % estimate circuit lifetime in days
        % cirlife_year = mttf_tddb(stbTc);
        % cirlife_day = cirlife_year * 365;
        cirlife_day = mttf_tddb(stbTc);
        
        % update total maintenance cost
        nodeC = c_bat / batlife_day + c_node / cirlife_day;
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

