% run experiments
clc;
clear;
close all;

% settings
run.IDSQ = false;
run.pSPIEL = true;
run.PSO = false;
run.ABC = false;
run.iter = 2;

% run experiments on small dataset
exp_small('pm2_5', run);
%exp_small('pm1', run);
%exp_small('pm10', run);
exp_small('humid', run);
%exp_small('temp', run);