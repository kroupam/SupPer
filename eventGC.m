function [value,isterminal,direction] = eventGC(~,y,par,parIn)
% This file is part of the program SupPer for the modelling of
% dynamic behavior and performance of Supercapacitors
% (Electrochemical Double Layer Capacitors)
% 
% Author: Martin Kroupa
% Oct 2015 - Mar 2016
% 
% Developed jointly at Imperial College London and University of Chemistry
% and Technology Prague

% Locate the time when the voltage passes through the boundaries

UMin = par.GC.UMin;
UMax = par.GC.UMax;

sigmaS = par.gen.sigmaS;

N = par.num.N;
Ns = par.num.Ns;
h = par.num.h;

js = parIn(1);

phi2End = y(2*N+Ns); % phi2 at the right current collector
phiDeltaEnd = y(4*N+Ns); % dphi at the right current collector
phi1End = phiDeltaEnd + phi2End; % phi1 at the right current collector
U = js*h/2.0/sigmaS + phi1End; % compute voltage


% and stop integration.
value(1) = UMax-U; % detect y-1/2 = 0
value(2) = UMin-U; % detect y-1/2 = 0
isterminal(1:2) = 1; % stop the integration
direction(1) = -1; % negative direction
direction(2) = 1; % negative direction