% set environmental variables to run ua:
uasource = '~/Ua/Ua2019/UaSource/';
addpath(genpath(uasource));
setenv('UaHomeDirectory',uasource);

% matlab utilities
addpath(genpath('/home/ronja/Ua/matlab_utilities/'), '-begin');

% PICO path
addpath(genpath('../'), '-begin');

