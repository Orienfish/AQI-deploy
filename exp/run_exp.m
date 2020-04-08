% run experiments
clc;
clear;
close all;

% settings
run.IDSQ = false;
run.pSPIEL = true;
run.PSO = true;
run.ABC = true;
run.iter = 1;

% run experiments on small dataset
exp_small('pm2_5', run);