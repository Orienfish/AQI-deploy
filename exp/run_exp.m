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
run.PSO = false;
run.ABC = false;
run.debugPlot = false;
run.iter = 10;

% run experiments on small dataset
exp_small('pm2_5', run, 10.0);
exp_small('pm1', run, 10.0);
exp_small('pm10', run, 10.0);
exp_small('humid', run, 10.0);
%exp_small('pm2_5', run, 7.0);
%exp_small('pm2_5', run, 8.0);
%exp_small('pm2_5', run, 9.0);
%exp_small('pm2_5', run, 11.0);
%exp_small('pm2_5', run, 12.0);
%exp_small('temp', run);

% run experiments on large dataset
%exp_large('pm2_5', run, 10.0);
%exp_large('pm1', run, 10.0);
%exp_large('pm10', run, 10.0);
%exp_large('humid', run, 10.0);