function [value,isterminal,direction] = eventCPC(~,y,par,parIn)
% This file is part of the program SupPer for the modelling of
% dynamic behavior and performance of Supercapacitors
% (Electrochemical Double Layer Capacitors)
% 
% Author: Martin Kroupa
% Oct 2015 - Mar 2016
% 
% Developed jointly at Imperial College London and University of Chemistry
% and Technology Prague

% Locate the time when height passes through zero in a decreasing
%direction

phiMin = par.CPC.phiMin;
phiMax = par.CPC.phiMax;

F = par.gen.F;
R = par.gen.R;
T = par.gen.T;
D = par.gen.D;
DS = par.gen.DS;
Ls = par.gen.Ls;

N = par.num.N;
h = par.num.h;

dim = par.spec.dim;

ps = parIn(1);

cN = y(2*N); %concentration at the interface between electrode and separator

switch dim
    case '1D'
        kappaL = F^2/R/T*2*(2*D*DS/(D+DS))*cN;
    case '1D1D'
        DMac = par.mes.DMac;
        kappaL = F^2/R/T*2*(2*DMac*DS/(DMac+DS))*cN;
end

% d_phi_fix = -jj*(h/2+wS/2)/kappaL; % fixed potetnial difference accros separator [V]
% U = 2*(y(N) - d_phi_fix);

U = y(N) + sqrt(y(N)^2 + ps*(Ls+h)/kappaL); % compute voltage

% and stop integration.
value(1) = phiMax-U; % detect y-1/2 = 0
value(2) = phiMin-U; % detect y-1/2 = 0
isterminal(1:2) = 1; % stop the integration
direction(1) = -1; % negative direction
direction(2) = 1; % negative direction