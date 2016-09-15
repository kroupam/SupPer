%% Parameters of Figures
function fig = parFigure()
% This file is part of the program SupPer for the modelling of
% dynamic behavior and performance of Supercapacitors
% (Electrochemical Double Layer Capacitors)
% 
% Author: Martin Kroupa
% Oct 2015 - Mar 2016
% 
% Developed jointly at Imperial College London and University of Chemistry
% and Technology Prague

% for working figures use the following settings:
fontSize = 14;
xSize = 17.0; % horizontal size of the figures
ySize = 12.5; % vertical size of the figures
line2 = 2; % thickness of lines

% % for presentations use the following settings:
% fontSize = 24;
% xSize = 22.2; %horizontal size of the figures
% ySize = 15.8; %vertical size of the figures
% line2 = 2; % thickness of lines

% Pack figure parameters
listFig = {'fieldNames' ...
    ,'fontSize' ...
    ,'xSize' ...
    ,'ySize' ...
    ,'line2'};

fig = v2struct(listFig);
end