% Test every functions
close all;
clear;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test battery lifetime function, plot out the lifetime at various temp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_bat = 10;
temp = linspace(-5, 40, n_bat);
batlife = zeros(1, n_bat);
for i = 1:n_bat
    batlife(i) = bat_lifetime(750, temp(i), 40, 0.1);
end
figure(1);
plot(temp, batlife);
title('Estimated battery time under various ambient temperature');
ylim([17 18.5]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test ambient temperature to core temperature conversion function
% plot out the core temperature at various ambient temperature
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_amb2core = 41;
Tamb = linspace(0, 40, n_amb2core);
Tcore = zeros(1, n_amb2core);
for i = 1:n_amb2core
    Tcore(i) = temp_amb2core(Tamb(i), 1.0);
end
figure(2);
plot(Tamb, Tcore);
title('Core temperature under various ambient temperature');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test MTTF_TDDB under various temperature and power
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_mttftemp = 41;
Tamb = linspace(0, 40, n_mttftemp);
n_pwr = 4;
pwr = linspace(1, 4, n_pwr);
MTTF = zeros(n_mttftemp, n_pwr);
for i = 1:n_mttftemp
    for j = 1:n_pwr
        MTTF(i, j) = mttf_tddb(1.1, Tamb(i), pwr(j));
    end
end
figure(3);
for j = 1:n_pwr
    plot(Tamb, MTTF(:, j));
    hold on;
end
hold off;
legend('pwr=1W', 'pwr=2W', 'pwr=3W', 'pwr=4W');
title('MTTF of TDDB under various ambient temperature');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test transmission power consumption under different distance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Pto = 0.52;       % 520mW
n_dist = 41;
dist = linspace(0, 10, n_dist);
Ptx = zeros(1, n_dist);
for i = 1:n_dist
    Ptx(i) = txPower(dist(i), Pto);
end
figure(4);
plot(dist, Ptx);
title('Average transmission power under various distance');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test average power consumption and tddb lifetime under different distance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Btx = 2500;       % 20kbps = 2500B/s
Brx = 2500;
Ltx = 1e3;        % 1kB
Lrx = 1e3;
Prx = 0.1;        % 100mW
Psen = 0.2;       % 200mW
tsen = 0.3;       % 300ms
Pslp = 52.2*1e-3; % 52.2mW
T = 10;           % 10s
avgPwr = zeros(1, n_dist);
MTTF_dist = zeros(1, n_dist);
for i = 1:n_dist
    avgPwr(i) = avgPower(dist(i), Btx, Ltx, Pto, Brx, Lrx, Prx, ...
        Psen, tsen, Pslp, T);
    MTTF_dist(i) = mttf_tddb(1.1, 25, avgPwr(i));
end
figure(5);
plot(dist, avgPwr);
title('Average total power under various distance');
figure(5);
plot(dist, MTTF_dist);
title('Average tddb lifetime under various distance');
