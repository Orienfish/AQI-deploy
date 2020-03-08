% Test every functions
% test battery lifetime functions, plot out the lifetime at various temp
temp = linspace(-5, 40, 10);
batlife = zeros(10);
for i = 1:10
    batlife(i) = bat_lifetime(750, temp(i), 40, 0.1);
end
figure();
plot(temp, batlife);
ylim([17 18.5])