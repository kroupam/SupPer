function dydt = rhsCV(t,y,par,parIn)
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
% Right Hand Side function for Cyclic Voltammetry simualtions
% t - time
% Y - vector of state variables
% DYDT - derivatives of state variables
% PARS - model parameters

%% Unpack parameters
UMin = par.CV.UMin;
UMax = par.CV.UMax;

dim = par.spec.dim;

nu = parIn(1);

t0 = (UMax - UMin)/nu; % half cycle period
U = UMax - nu*abs(mod(t,2*t0)-t0); % compute voltage

switch dim
    case '1D'
        dydt = model_1D(y,par,U);
    case '1D1D'
        dydt = model_1D1D(y,par,U);
end

end






















