function GC = parGC()
% This file is part of the program SupPer for the modelling of
% dynamic behavior and performance of Supercapacitors
% (Electrochemical Double Layer Capacitors)
% 
% Author: Martin Kroupa
% Oct 2015 - Mar 2016
% 
% Developed jointly at Imperial College London and University of Chemistry
% and Technology Prague

%% Parameters of Galvanostatic Cycling

NC = 20; % number of cycles
IVec = logspace(0,5,70); % current [A]
N_I = length(IVec); % length of the vector
UMin = 1.4; % [V] minimum voltage
UMax = 2.8; % [V] maximum voltage

% UMin = 0.01; % [V] minimum voltage
% UMax = 1.4 ; % [V] maximum voltage

%% Pack GC parameters
listGC = {'fieldNames' ...
    ,'NC' ...
    ,'IVec' ...
    ,'N_I' ...
    ,'UMin' ...
    ,'UMax' ...
    };

GC = v2struct(listGC);