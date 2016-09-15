function dydt = model_1D(y,par,U)
% This file is part of the program SupPer for the modelling of
% dynamic behavior and performance of Supercapacitors
% (Electrochemical Double Layer Capacitors)
% 
% Author: Martin Kroupa
% Oct 2015 - Mar 2016
% 
% Developed jointly at Imperial College London and University of Chemistry
% and Technology Prague

% 1D model of supercapacitor
% Y - vector of state variables (phi2 = potential in liquid phase ...
%                                dphi = potential difference (phi1-phi2) ...
%                                c = concentration of electrolyte ...
% DYDT - derivatives of state variables
% PAR - model parameters

% load parameters
R = par.gen.R;
F = par.gen.F;
epsTot = par.gen.epsTot;
aaTot = par.gen.aaTot;
% cB = par.gen.cB;
T = par.gen.T;
D = par.gen.D;
perm = par.gen.perm;
lambdaS = par.gen.lambdaS;
epsS = par.gen.epsS;
DS = par.gen.DS;
sigmaS = par.gen.sigmaS;

N = par.num.N;
h = par.num.h;
Ns = par.num.Ns;
hs = par.num.hs;

capMod = par.spec.capMod;

coeff = F^2/R/T*2*D; % coefficient for the computation of conductivity
coeffS = F^2/R/T*2*DS; % coefficient in separator

% read fields
phi2 = y(1:2*N+Ns); % potential in electrolyte (phi2, algebraic variable)
dphi = y(2*N+Ns+1:4*N+Ns); % potential difference (phi1-phi2)
c = y(4*N+Ns+1:6*N+2*Ns); % concentration of electrolyte

% extract phi1 (potential in the solid phase)
phi1 = NaN*phi2;
phi1(1:N) = dphi(1:N) + phi2(1:N);
phi1(N+Ns+1:2*N+Ns) = dphi(N+1:2*N) + phi2(N+Ns+1:2*N+Ns);

c(c<0.0) = 0.0; % concentration cannot be negative

% c(:) = cB; % supress temporarily the equation for concentration

%% calculate capacitance

phiDCap = phi1 - phi2; % the same as dphi, only with NaNs in separator (needed to satisfy format)

% calculate capacitance
switch capMod
    case 'constant'
        C = NaN*zeros(2*N+Ns,1);
        C(1:N) = perm/lambdaS;
        C(N+Ns+1:2*N+Ns) = perm/lambdaS;
    case 'GCS'
        lambdaD = sqrt(perm*R*T./(2*c*F^2)); % Debye length
        C_S = perm/lambdaS; % Stern layer capacitance
        C_GC = perm./lambdaD.*cosh(F*phiDCap/2.0/R/T); % potential-dependent capacitance
        C = 1.0./(1.0./C_GC + 1.0./C_S);
    case 'Bikerman'
        a = par.gen.a; % effective ion diameter [m]
        NA = par.gen.NA; % Avogadro number [mol^-1]
        nu = 2*a^3*c*NA; % packing fraction of ions
        lambdaD = sqrt(perm*R*T./(2*c*F^2)); % Debye length
        C = perm./lambdaD.*sinh(F*abs(phiDCap)/R/T)./( ...
            (1.0 + 2.0*nu.*sinh(F*phiDCap/2/R/T).^2).* ...
            sqrt(2.0./nu.*log(1.0 + 2.0*nu.*sinh(F*phiDCap/2/R/T).^2))...
            );
end

aC = aaTot*C;

%% calculate the time derivatives
dphi2dt = zeros(2*N+Ns,1); % initialize time derivative for phi2
ddphidt = zeros(2*N,1); % initialize time derivative for dphi
dcdt = zeros(2*N+Ns,1); % initialize time derivative for c

f2 = zeros(2*N+Ns+1,1); % initialize fluxes for phi2
f1 = zeros(2*N+Ns+1,1); % initialize fluxes for phi1
fc = zeros(2*N+Ns+1,1); % initialize fluxes for c

% left boundary (x = 0) (compute fluxes)
f2(1) = 0.0; % zero flux
f1(1) = -2.0*sigmaS/h*phi1(1); % zero potential at current collector
fc(1) = 0.0; % zero flux

% left electrode interior (compute fluxes)
for i=2:N 
    c_e = (c(i)+c(i-1))/2; % concentration at the edge
    kappaL = coeff*c_e;
    f2(i) = -kappaL/h*(phi2(i) - phi2(i-1)); % ionic current flux
    f1(i) = -sigmaS/h*(phi1(i) - phi1(i-1)); % ionic current flux
    fc(i) = -D/h*(c(i) - c(i-1)); % mass flux
end

% Boundary between electrode and separator (x = Le) (compute fluxes)
c_e = (h*c(N+1) + hs*c(N))/(h+hs); % concentration at the edge
kappaLS = coeffS*c_e; % on the side of separator
kappaL = coeff*c_e; % on the side of electrode
phi_e = (h*kappaLS*phi2(N+1) + hs*kappaL*phi2(N))/(hs*kappaL + h*kappaLS); % potential at the edge
f2(N+1) = -2.0*kappaL/h*(phi_e - phi2(N));

f1(N+1) = 0.0; % zero flux

c_e = (D*hs*c(N) + DS*h*c(N+1))/(D*hs + DS*h); % concentration at the edge
fc(N+1) = -2.0*D/h*(c_e - c(N)); % mass flux

% left electrode (compute balances)
for i=1:N
    ddphidt(i) = -1/aC(i)/h*(f2(i) - f2(i+1));
    dphi2dt(i) = 1/h*(f1(i) - f1(i+1)) + 1/h*(f2(i) - f2(i+1));
    dcdt(i) = 1.0/epsTot/h*(fc(i) - fc(i+1)) + aC(i)/2/epsTot/F*ddphidt(i);
end

% interior of separator (compute fluxes)
for i=N+2:N+Ns 
    c_e = (c(i)+c(i-1))/2; % concentration at the edge
    kappaL = coeffS*c_e;
    f2(i) = -kappaL/hs*(phi2(i) - phi2(i-1)); % ionic current flux
    fc(i) = -DS/hs*(c(i) - c(i-1)); % mass flux
end

% Boundary between electrode and separator (x = Le+Ls)
c_e = (h*c(N+Ns) + hs*c(N+Ns+1))/(h+hs); % concentration at the edge
kappaL = coeffS*c_e; % on the side of separator
kappaLE = coeff*c_e; % on the side of electrode
phi_e = (h*kappaL*phi2(N+Ns) + hs*kappaLE*phi2(N+Ns+1))/(hs*kappaLE + h*kappaL); % potential at the edge
f2(N+Ns+1) = -2.0*kappaL/hs*(phi_e - phi2(N+Ns));

f1(N+Ns+1) = 0.0; % zero flux

c_e = (D*hs*c(N+Ns+1) + DS*h*c(N+Ns))/(D*hs + DS*h); % concentration at the edge
fc(N+Ns+1) = -2.0*DS/hs*(c_e - c(N+Ns));

% separator (compute balances)
for i=N+1:N+Ns
    dphi2dt(i) = 1/hs*(f2(i) - f2(i+1));
    dcdt(i) = 1.0/epsS/hs*(fc(i) - fc(i+1));
end

% right electrode interior (compute fluxes)
for i=N+Ns+2:2*N+Ns 
    c_e = (c(i)+c(i-1))/2; % concentration at the edge
    kappaL = coeff*c_e;
    f2(i) = -kappaL/h*(phi2(i) - phi2(i-1)); % ionic current flux
    f1(i) = -sigmaS/h*(phi1(i) - phi1(i-1)); % ionic current flux
    fc(i) = -D/h*(c(i) - c(i-1)); % mass flux
end

% right boundary (x = 2*Le+Ls) (compute fluxes)
f2(2*N+Ns+1) = 0.0; % zero flux
f1(2*N+Ns+1) = -2.0*sigmaS/h*(U-phi1(2*N+Ns)); % total voltage at current collector
fc(2*N+Ns+1) = 0.0; % zero flux

% right electrode (compute balances)
for i=N+Ns+1:2*N+Ns
    ddphidt(i-Ns) = -1/aC(i)/h*(f2(i) - f2(i+1));
    dphi2dt(i) = 1/h*(f1(i) - f1(i+1)) + 1/h*(f2(i) - f2(i+1));
    dcdt(i) = 1.0/epsTot/h*(fc(i) - fc(i+1)) + aC(i)/2/epsTot/F*ddphidt(i-Ns);
end

% dcdt(:) = 0.0; % supress temporarily the equation for concentration

dydt = [dphi2dt; ddphidt; dcdt];

end

























