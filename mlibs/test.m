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
%n_mttftemp = 41;
%Tamb = linspace(0, 40, n_mttftemp);
%n_pwr = 4;
%pwr = linspace(1, 4, n_pwr);
%MTTF = zeros(n_mttftemp, n_pwr);
%for i = 1:n_mttftemp
%    for j = 1:n_pwr
%        MTTF(i, j) = mttf_tddb(1.1, Tamb(i), pwr(j));
%    end
%end
%figure(3);
%for j = 1:n_pwr
%    plot(Tamb, MTTF(:, j));
%    hold on;
%end
%hold off;
%legend('pwr=1W', 'pwr=2W', 'pwr=3W', 'pwr=4W');
%title('MTTF of TDDB under various ambient temperature');

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
% test average power consumption and tddb lifetime under different distance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
Tamb = linspace(0, 40, n_amb2core);
avgPwr = zeros(n_amb2core, n_dist);
%MTTF_dist = zeros(1, n_dist);
for k = 1:n_amb2core
    for i = 1:n_dist
        params.dtx = dist(i);
        [stbPwr, stbTc] = stbPower(params, Tamb(k));
        avgPwr(k, i) = stbPwr;
    %    MTTF_dist(i) = mttf_tddb(1.1, 25, avgPwr(i));
    end
end
figure(5);
for k = 1:n_amb2core
    plot(dist, avgPwr(k, :));
    hold on;
end
title('Average total power under various distance ambient temperature');
%figure(6);
%plot(dist, MTTF_dist);
%title('Average tddb lifetime under various distance');
