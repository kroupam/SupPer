function CV = parCV()
% This file is part of the program SupPer for the modelling of
% dynamic behavior and performance of Supercapacitors
% (Electrochemical Double Layer Capacitors)
% 
% Author: Martin Kroupa
% Oct 2015 - Mar 2016
% 
% Developed jointly at Imperial College London and University of Chemistry
% and Technology Prague

%% Parameters of Cyclic Voltammetry

NC = 20; % maximum number of cycles
nuVec = logspace(-4,3,70);
N_nu = length(nuVec); % length of the vector
UMin = 1.4; % [V] minimum voltage
UMax = 2.8; % [V] maximum voltage

% UMin = 0.01; % [V] minimum voltage
% UMax = 1.4 ; % [V] maximum voltage

%% Pack CV parameters
listCV = {'fieldNames' ...
    ,'NC' ...
    ,'nuVec' ...
    ,'N_nu' ...
    ,'UMin' ...
    ,'UMax' ...
    };

CV = v2struct(listCV);