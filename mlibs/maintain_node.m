function [nodeC, batlife_ratio, cirlife_ratio] = maintain_node(...
    params, curCenters, curProb)
%% Calculate the maintenance cost of a single node
%
% Args:
%   params: necessary parameters for power calculation
%   Centers: center of ambient temperature distribution
%   Prob: probability of ambient temperature in each bin around each center
%
% Return:
%   nodeC: the maintenance cost of this node
%   batlife_ratio: the battery lifetime in ratio of this node
%   cirlife_ratio: the circuit lifetime in ratio of this nods

nodeC = 0.0;
batlife_ratio = 0.0;
cirlife_ratio = 0.0;
% sum up through temperature distribution and
% calculate the expected maintenance cost of the current node
for j = 1:params.nbins
    % estimate power in W
    [stbPwr, stbTc] = stbPower(params, curCenters(j));

    % estimate battery lifetime in ratio
    % convert from W to mW then calculate average current draw
    I_mA = stbPwr * 1000 / params.Vdd;
    batlife_cur = bat_ratio(params.cap_bat, curCenters(j), I_mA, params.dt_bat_h);
    
    % estimate circuit lifetime in ratio
    cirlife_cur = mttf_ratio(stbTc);

    % update total maintenance cost
    nodeCur = params.c_bat / batlife_cur + params.c_node / cirlife_cur;

    % sum up expectations
    nodeC = nodeC + nodeCur * curProb(j);
    batlife_ratio = batlife_ratio + batlife_cur * curProb(j);
    cirlife_ratio = cirlife_ratio + cirlife_cur * curProb(j);
end      
    
end

