% Test every functions
% test battery lifetime function, plot out the lifetime at various temp
temp = linspace(-5, 40, 10);
batlife = zeros(10);
for i = 1:10
    batlife(i) = bat_lifetime(750, temp(i), 40, 0.1);
end
figure(1);
plot(temp, batlife);
ylim([17 18.5])

% test ambient temperature to core temperature conversion function
% plot out the core temperature at various ambient temperature
Tamb = linspace(20, 40, 20);
Tcore = zeros(20);
for i = 1:20
    Tcore(i) = temp_amb2core(Tamb(i), 1.0);
end
figure(2);
plot(Tamb, Tcore);