clear all; close all; clc  %#ok<*CLALL>
%% Program SupPer
% This is program SupPer for the modelling of dynamic behavior and
% performance of Supercapacitors (Electrochemical Double Layer Capacitors)
% 
% Author: Martin Kroupa
% Oct 2015 - Mar 2016
% 
% Developed jointly at Imperial College London and University of Chemistry
% and Technology Prague
%
%% Initialize (load default parameters)
fprintf('\n'); % new line in output
disp('********');
disp('Starting: SupPer''s Ready');
disp('********');

par.spec = parSpec; % load specifications about what should be computed
v2struct(par.spec); % unpack specification parameters
par = parameters(par);

%% Main computation

switch method
    case 'CV'
        par.CV = parCV; % load parameters of cyclic voltammetry
        res.CV = CV(par);
    case 'EIS'
        par.EIS = parEIS; % load parameters of electrochemical impedance spectroscopy
        res.EIS = EIS(par);
    case 'GC'
        par.GC = parGC; % load parameters of galvanostatic cycling
        res.GC = GC(par);
    otherwise
        error(['SupPer: ' method ' is an unknown method'])
end

%% Plot current and voltage characteristics

k = 60; % index of frequency / scan rate to plot

par.fig = parFigure; % load parameters for figures
figureCurVolTime(res.(method), par, k)

%% Plot capacity vs. frequency

par.fig = parFigure; % load parameters for figures
figureCapFreq(res, par)

%% Plot spatial distributions

k = 10; % index of frequency / scan rate to plot
tBeg = 0; % initial time
tEnd = 1000.0; % end time

% plot results
figurePotConXTime(res.(method), par, k, tBeg, tEnd)

switch dim
    case '1D1D'
        figurePotConXY(res.(method), par, k, tEnd)
end

disp('Finished: All graphs plotted, SupPer is over');









