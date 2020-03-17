function [avgPwr] = avgPower(dtx, Btx, Ltx, Brx, Lrx, Prx, Pslp, T)
% Calculate estimated power for a sensor.
%
% Args:
% (i) for transmission: distance in m(dtx), bandwidth in B/s (Btx), data size in B (Ltx)
% (ii) for receiving: bandwidth in B/s (Brx), data size in B (Lrx), avg. recv. power in W (Prx)
% (iii) for sleep: avg. sleep power in W (Pslp)
% (iv) time for one period or data sampling period in s (T)
%
% Return:
%   avgPwr: average power consumption in W

% calculate avg. transmission power and time
Ptx = txPower(dtx);
ttx = Ltx / Btx;
% fprintf('tx power: %f tx time: %f\n', Ptx, ttx);

% calculate avg. receiving time
trx = Lrx / Brx;
% fprintf('rx power: %f rx time: %f\n', Prx, trx);

% calculate total energy and power in one period
energy = Ptx * ttx + Prx * trx + Pslp * (T - ttx - trx);
% fprintf('tx energy: %f rx energy: %f sleep energy: %f\n', ...
%    Ptx * ttx, Prx * trx, Pslp * (T - ttx - trx));
avgPwr = energy / T;
end

