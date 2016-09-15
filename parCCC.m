function CCC = parCCC()
% This file is part of the program SupPer for the modelling of
% dynamic behavior and performance of Supercapacitors
% (Electrochemical Double Layer Capacitors)
% 
% Author: Martin Kroupa
% Oct 2015 - Mar 2016
% 
% Developed jointly at Imperial College London and University of Chemistry
% and Technology Prague

%% Parameters of Constant Current Charging followed by Constant Voltage Discharging

I_fix = 100; % constant charging current [A]
UMin = 1.6; % [V] minimum voltage
UDis = 1.4; % [V] discharge voltage

%% Pack GC parameters
listCCC = {'fieldNames' ...
    ,'I_fix' ...
    ,'UMin' ...
    ,'UDis' ...
    };

CCC = v2struct(listCCC);