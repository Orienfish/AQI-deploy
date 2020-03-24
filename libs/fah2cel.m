function [T_cel] = fah2cel(T_fah)
%% Convert Fahrenheit to Celsius.
% Args:
%   T_fah: array of temperature in Fahrenheit
%
% Return:
%   T_cel: array of temperature in Celsius

T_cel = (5/9)*(T_fah - 32);
end

