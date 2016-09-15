function dydt = rhsEIS(t,y,par,parIn)
% This file is part of the program SupPer for the modelling of
% dynamic behavior and performance of Supercapacitors
% (Electrochemical Double Layer Capacitors)
% 
% Author: Martin Kroupa
% Oct 2015 - Mar 2016
% 
% Developed jointly at Imperial College London and University of Chemistry
% and Technology Prague
% 
% Right Hand Side function for Electrochemical Impedance Spectroscopy
% t - time
% Y - vector of state variables
% DYDT - derivatives of state variables
% PAR - model parameters

%% Unpack parameters
UAmp = par.EIS.UAmp;
UFix = par.EIS.UFix;

dim = par.spec.dim;

f = parIn(1);

U = UFix + UAmp*exp(1i*2*pi*f*t); % compute voltage

switch dim
    case '1D'
        dydt = model_1D(y,par,U);
    case '1D1D'
        dydt = model_1D1D(y,par,U);
end

end






















