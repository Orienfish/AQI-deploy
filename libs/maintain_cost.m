function [C] = maintain_cost(Xa, Ta, commMST)
%% Compute the maintenance cost from a given deployment.
%
% Args:
%   Xa: a given deployment list, [lat lon]
%   Ta: average temperature estimation at Xa in Celsius 
%   commMST: generated feasible MST for connnection
%
% Return:
%   C: the maintenance cost

% settings for sensor workloads
Btx = 2500;    % 20kbps = 2500B/s
Brx = 2500;
Ltx = 1e3;     % 1kB
Lrx = 1e3;
Prx = 22.2e-3; % 22.2mW
Pslp = 1e-4;   % 0.1mW
T = 10;        % sampling interval in seconds

% settings for battery
cap_bat = 750; % initial battery capacity in mAh
dt_bat_h = 1;  % time resolution of battery in hours
V = 3.3;       % supply voltage
c_bat = 1;    % cost to replace battery

% setting for circuit
V_gs = 1.1;    % gate voltage
c_node = 100;  % cost to replace node
C = 0;         % total maintenance cost

% iterative through every node in Xa
for i = 1:size(Xa, 1)
    % estimate power in W
    txDist_km = commMST(~isnan(commMST(:, i)), i); % get the unique non-nan
    txDist_m = txDist_km / 1000;  % convert to meters
    if length(txDist_m) ~= 1
        error('Invalid commMST!'); % uniqueness check
    end
    child_cnt = nnz(~isnan(commMST(i, :))); % get child cnt of node i
    avgPwr = avgPower(txDist_m, Btx, Ltx, Brx, child_cnt * Lrx, Prx, Pslp, T);
    fprintf('avg pwr of node [%f %f] is %f W\n', Xa(i, 1), Xa(i, 2), avgPwr);
    
    % estimate battery lifetime in days
    I_mA = avgPwr * 1000 / V; % convert from W to mW then calculate average current draw
    batlife_h = bat_lifetime(cap_bat, Ta(i), I_mA, dt_bat_h);
    batlife_day = batlife_h / 24;
    fprintf('battery life of node [%f %f] is %f days\n', Xa(i, 1), Xa(i, 2), batlife_day);
    
    % estimate circuit lifetime in days
    cirlife_year = mttf_tddb(V_gs, Ta(i), avgPwr);
    cirlife_day = cirlife_year * 365;
    fprintf('circuit life of node [%f %f] is %f days\n', Xa(i, 1), Xa(i, 2), cirlife_day);
    
    % update total maintenance cost
    nodeC = c_bat / batlife_day + c_node / cirlife_day;
    C = C + nodeC;
    fprintf('maintenance cost of node [%f %f] is %f\n', Xa(i, 1), Xa(i, 2), nodeC);
end
end

