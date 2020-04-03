% Test every functions
clc;
clear;
close all;
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
n_amb2core = 5;
iter = 400;
Tamb = linspace(0, 40, n_amb2core);
Tcore = zeros(n_amb2core, iter);
for i = 1:n_amb2core
    Tcore(i, 1) = Tamb(i);
    for j = 2:iter
        Tcore(i, j) = temp_amb2core(Tamb(i), 1.0, Tcore(i, j-1));
    end
end
figure(2);
for i = 1:n_amb2core
    plot(Tcore(i, :));
    hold on;
end
title('Core temperature under various ambient temperature');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test MTTF_TDDB under various temperature and power
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_mttftemp = 41;
Tc = linspace(0, 100, n_mttftemp); % core temperature in Celsius
n_pwr = 4;
MTTF = zeros(1, n_mttftemp);
for i = 1:n_mttftemp
    MTTF(i) = mttf_tddb(Tc(i) + 273.15); % convert to Kelvin
end
figure(3);
plot(Tc, MTTF);
title('MTTF of TDDB under various core temperature');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test transmission power consumption under different distance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params.Pto = 0.32;       % 320mW
n_dist = 41;
dist = linspace(0, 10, n_dist);
Ptx = zeros(1, n_dist);
for i = 1:n_dist
    Ptx(i) = txPower(dist(i), params.Pto);
end
figure(4);
plot(dist, Ptx);
title('Average transmission power under various distance');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test average power consumption, battery lifetime and tddb lifetime 
% under different distance and ambient temperature
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_amb = 5;
Tamb = linspace(0, 40, n_amb);
n_dist = 41;
dist = linspace(0, 10, n_dist);
avgPwr = zeros(n_amb, n_dist);
MTTF_dist = zeros(n_amb, n_dist);
batlife_dist = zeros(n_amb, n_dist);
% settings for power
params.Btx = 2500;       % 20kbps = 2500B/s
params.Brx = 2500;
params.Ltx = 1e3;        % 1kB
params.Lrx = 1e3;
params.Prx = 0.1;        % 100mW
params.Psen = 0.2;       % 200mW
params.tsen = 0.3;       % 300ms
params.T = 10;           % 10s
params.Vdd = 3.3;          % 3.3v
params.f = 300e6;        % 300MHz
% settings for battery
cap_bat = 20000;   % initial battery capacity in mAh
dt_bat_h = 1;     % time resolution of battery in hours

for k = 1:n_amb
    for i = 1:n_dist
        params.dtx = dist(i);
        [stbPwr, stbTc] = stbPower(params, Tamb(k));
        
        avgPwr(k, i) = stbPwr;
        MTTF_dist(k, i) = mttf_tddb(stbTc);
        
        % convert from W to mW then calculate average current draw
        I_mA = stbPwr * 1000 / params.Vdd; 
        batlife_h = bat_lifetime(cap_bat, Tamb(k), I_mA, dt_bat_h);
        batlife_dist(k, i) = batlife_h;
    end
end
figure(5);
for k = 1:n_amb
    plot(dist, avgPwr(k, :));
    hold on;
end
title('Average total power under various distance and ambient temperature');
figure(6);
for k = 1:n_amb
    plot(dist, MTTF_dist(k, :));
    hold on;
end
title('Average tddb lifetime under various distance and ambient temperature');
figure(7);
for k = 1:n_amb
    plot(dist, batlife_dist(k, :));
    hold on;
end
title('Average battery lifetime under various distance and ambient temperature');
