function dydt = rhsGC(~,y,par,parIn)
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
% Right Hand Side function for Galvanostatic Cycling simualtions
% t - time
% Y - vector of state variables
% DYDT - derivatives of state variables
% PAR - model parameters

%% Unpack parameters
sigmaS = par.gen.sigmaS;

N = par.num.N;
Ns = par.num.Ns;
h = par.num.h;

dim = par.spec.dim;

js = parIn(1);

phi2End = y(2*N+Ns); % phi2 at the right current collector
phiDeltaEnd = y(4*N+Ns); % dphi at the right current collector
phi1End = phiDeltaEnd + phi2End; % phi1 at the right current collector
U = js*h/2.0/sigmaS + phi1End; % compute voltage

switch dim
    case '1D'
        dydt = model_1D(y,par,U);
    case '1D1D'
        dydt = model_1D1D(y,par,U);
end

end






















