% A small matlab wrapper script for simulating the controlled unicycle
% trajectory
clear
close all

% Load information (numAbs, numCont)
run('C/controller_info.m')

% Plot state space
unicycle('S', numAbs, numCont, 'f');
close all;

% Plot obstacles and goal
unicycle('P', numAbs, numCont, 'f');
close all;

% Simulate controlled trajectory
unicycle('R', numAbs, numCont, 'f');
