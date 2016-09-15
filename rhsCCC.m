function dydt = rhsCCC(~,y,par,parIn)
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
% Right Hand Side function for Constant Current Charging
% t - time
% Y - vector of state variables
% DYDT - derivatives of state variables
% PAR - model parameters

%% Unpack parameters
sigmaS = par.gen.sigmaS;

N = par.num.N;
Ns = par.num.Ns;
h = par.num.h;

js = parIn(1);

phi2End = y(2*N+Ns);
phiDeltaEnd = y(4*N+Ns);
phi1End = phiDeltaEnd + phi2End;
U = js*h/2.0/sigmaS + phi1End;

% kappaL = F^2/R/T*2*DS*cB;
% d_phi_fix = -js*hs/2.0/kappaL; % fixed potetnial difference accros separator [V]
% 
% U = 2*(y(N+Ns) - d_phi_fix);

dydt = model_1D(y,par,U);

end






















