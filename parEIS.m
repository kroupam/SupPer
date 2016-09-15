function EIS = parEIS()
% This file is part of the program SupPer for the modelling of
% dynamic behavior and performance of Supercapacitors
% (Electrochemical Double Layer Capacitors)
% 
% Author: Martin Kroupa
% Oct 2015 - Mar 2016
% 
% Developed jointly at Imperial College London and University of Chemistry
% and Technology Prague

%% Parameters of Electrochemical Impedance Spectroscopy

fVec = logspace(-4,3,70); % frequency [1/s]
N_f = length(fVec); % length of the vector
UAmp = 0.005; % [V] voltage amplitude
UFix = 2.1; % [V] fixed voltage

nT = 10; % number of periods to be simulated

%% Pack EIS parameters
listEIS = {'fieldNames' ...
    ,'fVec' ...
    ,'N_f' ...
    ,'UAmp' ...
    ,'UFix' ...
    ,'nT' ...
    };

EIS = v2struct(listEIS);