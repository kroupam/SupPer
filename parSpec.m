function spec = parSpec()
% This file is part of the program SupPer for the modelling of
% dynamic behavior and performance of Supercapacitors
% (Electrochemical Double Layer Capacitors)
% 
% Author: Martin Kroupa
% Oct 2015 - Mar 2016
% 
% Developed jointly at Imperial College London and University of Chemistry
% and Technology Prague

% specification of the simulation
%% Choose Analytical Method
method = 'CV';
% method = 'EIS';
% method = 'GC';

%% Choose Dimension
dim = '1D';
% dim = '1D1D';

%% Choose Model for Capacitance
capMod = 'constant';
% capMod = 'GCS';
% capMod = 'Bikerman';

%% Test acceptability and print information

disp('Chosen Method:')
switch method
    case 'CV'
        disp('Cyclic Voltammetry')
    case 'EIS'
        disp('Electrochemical Impedance Spectroscopy')
    case 'GC'
        disp('Galvanostatic Cycling')
    otherwise
        error(['SupPer: ' method ' is an unknown method'])
end
fprintf('\n'); % new line in output

disp('Chosen Dimension:')
switch dim
    case '1D'
        disp('1D')
    case '1D1D'
        disp('1D + 1D (Warning: the simulation will take longer to compute)')
    otherwise
        error(['SupPer: ' dim ' is an unknown dimension'])
end
fprintf('\n'); % new line in output

disp('Chosen Model for Capacitance:')
switch capMod
    case 'constant'
        disp('constant')
    case 'GCS'
        disp('Gouy-Chapman-Stern')
    case 'Bikerman'
        disp('Bikerman')
    otherwise
        error(['SupPer: ' capMod ' is an unknown model'])
end
fprintf('\n'); % new line in output

%% Pack specification parameters
listSpec = {'fieldNames' ...
    ,'method' ...
    ,'dim' ...
    ,'capMod' ...
    };

spec = v2struct(listSpec);


