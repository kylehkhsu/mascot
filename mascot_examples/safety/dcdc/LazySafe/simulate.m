% A small matlab wrapper script for simulating the controlled unicycle
% trajectory
clear
close all

% Load information (numAbs, numCont)
run('C/controller_info.m')

% Plot state space
dcdc('S', numAbs, numCont, 'f');
close all;

% % Plot obstacles and goal
% dcdc('P', numAbs, numCont, 'f');
% close all;

% Simulate controlled trajectory
dcdc('Safe', numAbs, numCont, 'f');
