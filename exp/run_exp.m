% run experiments
clc;
clear;
close all;
addpath('../libs/');
addpath('../mlibs/');
addpath('../lldistkm/');
addpath('../gp/');
addpath('../SFO/sfo/');

% settings
run.IDSQ = false;
run.pSPIEL = false;
run.PSO = true;
run.ABC = false;
run.debugPlot = false;
run.iter = 10;

% run experiments on small dataset
%exp_small('pm2_5', run, 10.0, 16);
%exp_small('pm1', run, 10.0, 16);
%exp_small('pm10', run, 10.0, 16);
%exp_small('humid', run, 10.0, 16);
%exp_small('temp', run, 10.0, 16); % cannot run pSPIEL, cov_vd matrix is not PSD
%exp_small('pm2_5', run, 7.0, 16);
%exp_small('pm2_5', run, 8.0, 16);
%exp_small('pm2_5', run, 9.0, 16);
%exp_small('pm2_5', run, 11.0, 16);
%exp_small('pm2_5', run, 12.0, 16);

% test execution time with various number of sensors
%exp_small('pm2_5', run, 10.0, 7);
%exp_small('pm2_5', run, 10.0, 9);
%exp_small('pm2_5', run, 10.0, 11);
%exp_small('pm2_5', run, 10.0, 13);
%exp_small('pm2_5', run, 10.0, 15);


% run experiments on large dataset
% Quota -> # of sensors in pSPIEL
% 50 -> 85, 40 -> 62
exp_large('pm2_5', run, 40.0, 54); 
exp_large('pm1', run, 40.0, 54);
exp_large('pm10', run, 40.0, 54);
exp_large('humid', run, 40.0, 54);
exp_large('temp', run, 40.0, 54);
exp_large('pm2_5', run, 32.0, 54);
exp_large('pm2_5', run, 34.0, 54);
exp_large('pm2_5', run, 36.0, 54);
exp_large('pm2_5', run, 38.0, 54);