function [C] = maintain_cost(Xa, Ta, commMST, logging)
%% Compute the maintenance cost from a given deployment.
%
% Args:
%   Xa: a given deployment list, [lat lon]
%   Ta: average temperature estimation at Xa in Celsius 
%   commMST: generated feasible MST for connnection
%   logging: logging flag
%
% Return:
%   C: the maintenance cost

% settings for sensor workloads
Pto = 0.52;    % 520mW
Btx = 2500;       % 20kbps = 2500B/s
Brx = 2500;
Ltx = 1e3;        % 1kB
Lrx = 1e3;
Prx = 0.1;        % 100mW
Psen = 0.2;       % 200mW
tsen = 0.3;       % 300ms
Pslp = 52.2*1e-3; % 52.2mW
T = 10;           % 10s

% settings for battery
cap_bat = 2000;   % initial battery capacity in mAh
dt_bat_h = 1;     % time resolution of battery in hours
V = 3.3;          % supply voltage
c_bat = 1;        % cost to replace battery

% setting for circuit
V_gs = 1.1;       % gate voltage
c_node = 100;     % cost to replace node
C = 0;            % total maintenance cost

% iterative through every node in Xa
for i = 1:size(Xa, 1)
    % estimate power in W
    txDist_km = commMST(~isnan(commMST(:, i)), i); % get the unique non-nan
    txDist_m = txDist_km / 1000;  % convert to meters
    if length(txDist_m) ~= 1
        error('Invalid commMST!'); % uniqueness check
    end
    child_cnt = nnz(~isnan(commMST(i, :))); % get child cnt of node i
    avgPwr = avgPower(txDist_m, Btx, Ltx, Pto, Brx, child_cnt * Lrx, Prx, ...
        Psen, tsen, Pslp, T);
    
    % estimate battery lifetime in days
    I_mA = avgPwr * 1000 / V; % convert from W to mW then calculate average current draw
    batlife_h = bat_lifetime(cap_bat, Ta(i), I_mA, dt_bat_h);
    batlife_day = batlife_h / 24;
    
    % estimate circuit lifetime in days
    cirlife_year = mttf_tddb(V_gs, Ta(i), avgPwr);
    cirlife_day = cirlife_year * 365;
    
    % update total maintenance cost
    nodeC = c_bat / batlife_day + c_node / cirlife_day;
    C = C + nodeC;
    
    if logging
        disp(['node [' num2str(Xa(i, 1)) ' ' num2str(Xa(i, 2))]);
        disp(['amb temp: ' num2str(Ta(i)) 'avg pwr: ' num2str(avgPwr)]);
        disp(['bat life: ' batlife_day 'cir life: ' cirlife_day 'main cost: ' nodeC]);
    end
end
end

