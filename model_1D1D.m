function dydt = model_1D1D(y,par,U)
% This file is part of the program SupPer for the modelling of
% dynamic behavior and performance of Supercapacitors
% (Electrochemical Double Layer Capacitors)
% 
% Author: Martin Kroupa
% Oct 2015 - Mar 2016
% 
% Developed jointly at Imperial College London and University of Chemistry
% and Technology Prague

% 1D+1D model of supercapacitor
% Y - vector of state variables (phi2 = potential in liquid phase in macro-pores ...
%                                dphi = potential difference (phi1-phi2) ...
%                                c = concentration of electrolyte in macro-pores ...
%                                dphi_p = potential difference (phi1-phi_p) in meso-pores ...
%                                cp = concentration of electrolyte in meso-pores)
% DYDT - derivatives of state variables
% PAR - model parameters

% load parameters
R = par.gen.R;
F = par.gen.F;
% cB = par.gen.cB;
T = par.gen.T;
perm = par.gen.perm;
lambdaS = par.gen.lambdaS;
epsS = par.gen.epsS;
DS = par.gen.DS;
sigmaS = par.gen.sigmaS;

N = par.num.N;
h = par.num.h;
Ns = par.num.Ns;
hs = par.num.hs;
Np = par.num.Np;
hp = par.num.hp;
Sr = par.num.Sr;
V = par.num.V;

epsMac = par.mes.epsMac;
epsMes = par.mes.epsMes;
DMac = par.mes.DMac;
DMes = par.mes.DMes;
% Rp = par.mes.Rp;
aaMac = par.mes.aaMac;
aaMes = par.mes.aaMes;

capMod = par.spec.capMod;

coeffMac = F^2/R/T*2*DMac; % coefficient for the computation of conductivity (macro-pores)
coeffMes = F^2/R/T*2*DMes; % coefficient for the computation of conductivity (meso-pores)
coeffS = F^2/R/T*2*DS; % coefficient in separator

% read fields
phi2 = y(1:2*N+Ns); % potential in electrolyte (phi2, algebraic variable)
dphi = y(2*N+Ns+1:4*N+Ns); % potential difference (phi1-phi2)
c = y(4*N+Ns+1:6*N+2*Ns); % concentration of electrolyte in macro-pores
dphi_p = y(6*N+2*Ns + 1:6*N+2*Ns + 2*N*Np); % potential difference (phi1-phi2p) in meso-pores
cp = y(6*N+2*Ns + 2*N*Np + 1:6*N+2*Ns + 4*N*Np); % concentration of electrolyte in meso-pores

% extract phi1 (potential in the solid phase)
phi1 = NaN*phi2;
phi1(1:N) = dphi(1:N) + phi2(1:N);
phi1(N+Ns+1:2*N+Ns) = dphi(N+1:2*N) + phi2(N+Ns+1:2*N+Ns);

% extract phi2_p (potential in the liquid phase in meso-pores)
phi2_p = NaN*dphi_p;
for i=1:N
    phi2_p((i-1)*Np+1:(i-1)*Np+Np) = phi1(i) - dphi_p((i-1)*Np+1:(i-1)*Np+Np);
end
for i=N+1:2*N
    phi2_p((i-1)*Np+1:(i-1)*Np+Np) = phi1(Ns+i) - dphi_p((i-1)*Np+1:(i-1)*Np+Np);
end

c(c<0.0) = 0.0; % concentration cannot be negative
cp(cp<0.0) = 0.0; % concentration cannot be negative

% % supress temporarily the equations for potential
% phi(:) = 0.0; 
% phi_p(:) = 0.0;
% U = 0.0;

% % supress temporarily the equations for concentration
% c(:) = cB; 
% cp(:) = cB;

%% calculate the macro-pore phase
dphi2dt = zeros(2*N+Ns,1); % initialize time derivative for phi2
ddphidt = zeros(2*N,1); % initialize time derivative for dphi
dcdt = zeros(2*N+Ns,1); % initialize time derivative for c

phiDCap = phi1 - phi2; % the same as dphi, only with NaNs in separator (needed to satisfy format)

% calculate capacitance
switch capMod
    case 'constant'
        CMac = NaN*zeros(2*N+Ns,1);
        CMac(1:N) = perm/lambdaS;
        CMac(N+Ns+1:2*N+Ns) = perm/lambdaS;
    case 'GCS'
        lambdaD = sqrt(perm*R*T./(2*c(1:N)*F^2)); % Debye length
        C_S = perm/lambdaS; % Stern layer capacitance
        C_GC = perm./lambdaD.*cosh(F*phiDCap/2.0/R/T); % potential-dependent capacitance
        CMac = 1.0./(1.0./C_GC + 1.0./C_S);
    case 'Bikerman'
        a = par.gen.a; % effective ion diameter [m]
        NA = par.gen.NA; % Avogadro number [mol^-1]
        nu = 2*a^3*c(1:N)*NA; % packing fraction of ions
        lambdaD = sqrt(perm*R*T./(2*c(1:N)*F^2)); % Debye length
        CMac = perm./lambdaD.*sinh(F*abs(phiDCap)/R/T)./( ...
            (1.0 + 2.0*nu.*sinh(F*phiDCap/2/R/T).^2).* ...
            sqrt(2.0./nu.*log(1.0 + 2.0*nu.*sinh(F*phiDCap/2/R/T).^2))...
            );
end

aCMac = aaMac*CMac;

f2 = zeros(2*N+Ns+1,1); % initialize fluxes for phi2
f1 = zeros(2*N+Ns+1,1); % initialize fluxes for phi1
fc = zeros(2*N+Ns+1,1);  % initialize fluxes for c

% left boundary (x = 0) (compute fluxes)
f2(1) = 0.0; % zero flux
f1(1) = -2.0*sigmaS/h*phi1(1); % zero potential at current collector
fc(1) = 0.0; % zero flux

% left electrode interior (compute fluxes)
for i=2:N 
    c_e = (c(i)+c(i-1))/2; % concentration at the edge
    kappaL = coeffMac*c_e;
    f2(i) = -kappaL/h*(phi2(i) - phi2(i-1)); % ionic current flux
    f1(i) = -sigmaS/h*(phi1(i) - phi1(i-1)); % ionic current flux
    fc(i) = -DMes/h*(c(i) - c(i-1)); % mass flux
end

% Boundary between electrode and separator (x = Le) (compute fluxes)
c_e = (h*c(N+1) + hs*c(N))/(h+hs); % concentration at the edge
kappaLS = coeffS*c_e; % on the side of separator
kappaL = coeffMac*c_e; % on the side of electrode
phi_e = (h*kappaLS*phi2(N+1) + hs*kappaL*phi2(N))/(hs*kappaL + h*kappaLS); % potential at the edge
f2(N+1) = -2.0*kappaL/h*(phi_e - phi2(N));

f1(N+1) = 0.0; % zero flux

c_e = (DMac*hs*c(N) + DS*h*c(N+1))/(DMac*hs + DS*h); % concentration at the edge
fc(N+1) = -2.0*DMac/h*(c_e - c(N)); % mass flux

% left electrode (compute balances)
for i=1:N
    kappaL = coeffMes*c(i);
    
    ddphidt(i) = -1/aCMac(i)/h*(f2(i) - f2(i+1)) + ...
        aaMac/aCMac(i)*(2*kappaL*(phi2(i) - phi2_p((i-1)*Np + Np))/hp);
    
    dphi2dt(i) = 1/h*(f1(i) - f1(i+1)) + 1/h*(f2(i) - f2(i+1));
    
    dcdt(i) = 1.0/epsMac/h*(fc(i) - fc(i+1)) + aCMac(i)/2/epsMac/F*ddphidt(i) - ...
        aaMac/epsMac*(2*DMes*(c(i) - cp((i-1)*Np + Np))/hp);
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
kappaLE = coeffMac*c_e; % on the side of electrode
phi_e = (h*kappaL*phi2(N+Ns) + hs*kappaLE*phi2(N+Ns+1))/(hs*kappaLE + h*kappaL); % potential at the edge
f2(N+Ns+1) = -2.0*kappaL/hs*(phi_e - phi2(N+Ns));

f1(N+Ns+1) = 0.0; % zero flux

c_e = (DMac*hs*c(N+Ns+1) + DS*h*c(N+Ns))/(DMac*hs + DS*h); % concentration at the edge
fc(N+Ns+1) = -2.0*DS/hs*(c_e - c(N+Ns));

% separator (compute balances)
for i=N+1:N+Ns
    dphi2dt(i) = 1/hs*(f2(i) - f2(i+1));
    dcdt(i) = 1.0/epsS/hs*(fc(i) - fc(i+1));
end

% right electrode interior (compute fluxes)
for i=N+Ns+2:2*N+Ns
    c_e = (c(i)+c(i-1))/2; % concentration at the edge
    kappaL = coeffMac*c_e;
    f2(i) = -kappaL/h*(phi2(i) - phi2(i-1)); % ionic current flux
    f1(i) = -sigmaS/h*(phi1(i) - phi1(i-1)); % ionic current flux
    fc(i) = -DMac/h*(c(i) - c(i-1)); % mass flux
end

% right boundary (x = 2*Le+Ls) (compute fluxes)
f2(2*N+Ns+1) = 0.0; % zero flux
f1(2*N+Ns+1) = -2.0*sigmaS/h*(U-phi1(2*N+Ns)); % total voltage at current collector
fc(2*N+Ns+1) = 0.0; % zero flux

% right electrode (compute balances)
for i=N+Ns+1:2*N+Ns
    kappaL = coeffMes*c(i);
    
    ddphidt(i-Ns) = -1/aCMac(i)/h*(f2(i) - f2(i+1)) + ...
        aaMac/aCMac(i)*(2*kappaL*(phi2(i) - phi2_p((i-1-Ns)*Np + Np))/hp);
    
    dphi2dt(i) = 1/h*(f1(i) - f1(i+1)) + 1/h*(f2(i) - f2(i+1));
    
    dcdt(i) = 1.0/epsMac/h*(fc(i) - fc(i+1)) + aCMac(i)/2/epsMac/F*ddphidt(i-Ns) - ...
        aaMac/epsMac*(2*DMes*(c(i) - cp((i-1-Ns)*Np + Np))/hp);
end

%% calculate the meso-pore phase
ddphi_pdt = zeros(2*N*Np,1); % initialize time derivative for dphi_p
dcpdt = zeros(2*N*Np,1); % initialize time derivative for cp

for i=1:2*N
    
    if i>=N+1
        iAdd = Ns; % add to satisfy format
    else
        iAdd = 0;
    end
    
    fp = zeros(Np+1,1); % initialize fluxes for phi_p
    fc = zeros(Np+1,1); % initialize fluxes for cp
    
    phiDCap = phi1(i+iAdd) - phi2_p((i-1)*Np+1 : (i-1)*Np+Np); % the same as dphi_p, only with NaNs in separator (needed to satisfy format)
    
    % calculate capacitance
    switch capMod
        case 'constant'
            CMes = zeros(Np,1);
            CMes(:) = perm/lambdaS;
        case 'GCS'
            lambdaD = sqrt(perm*R*T./(2*cp((i-1)*Np+1 : (i-1)*Np+Np)*F^2)); % Debye length
            C_S = perm/lambdaS; % Stern layer capacitance
            C_GC = perm./lambdaD.*cosh(F*phiDCap/2.0/R/T); % potential-dependent capacitance
            CMes = 1.0./(1.0./C_GC + 1.0./C_S);
        case 'Bikerman'
            nu = 2*a^3*cp((i-1)*Np+1 : (i-1)*Np+Np)*NA; % packing fraction of ions
            lambdaD = sqrt(perm*R*T./(2*cp((i-1)*Np+1 : (i-1)*Np+Np)*F^2)); % Debye length
            CMes = perm./lambdaD.*sinh(F*abs(phiDCap)/R/T)./( ...
                (1.0 + 2.0*nu.*sinh(F*phiDCap/2/R/T).^2).* ...
                sqrt(2.0./nu.*log(1.0 + 2.0*nu.*sinh(F*phiDCap/2/R/T).^2))...
                );
    end
    
    aCMes = aaMes*CMes;
    
    % centre of the particle (r = 0) (compute fluxes)
    fp(1) = 0.0; % zero flux in the centre
    fc(1) = 0.0; % zero flux in the centre
    
    c_e = (0.5*cp((i-1)*Np + 1) + cp((i-1)*Np + 2))/1.5; % concentration at the edge
    kappaL = coeffMes*c_e;
    fp(2) = -kappaL/(1.5*hp)*(phi2_p((i-1)*Np + 2) - phi2_p((i-1)*Np + 1)); % current flux
    fc(2) = -DMes/(1.5*hp)*(cp((i-1)*Np + 2) - cp((i-1)*Np + 1)); % mass flux
    
    % interior of the particle (compute fluxes)
    for j=3:Np
        c_e = (cp((i-1)*Np + j-1) + cp((i-1)*Np + j))/2; % concentration at the edge
        kappaL = coeffMes*c_e;
        fp(j) = -kappaL/hp*(phi2_p((i-1)*Np + j) - phi2_p((i-1)*Np + j-1)); % current flux
        fc(j) = -DMes/hp*(cp((i-1)*Np + j) - cp((i-1)*Np + j-1)); % mass flux
    end
    
    % surface of the particle (r = Rp) (compute fluxes)
    c_e = c(i+iAdd); % concentration at the edge
    kappaL = coeffMes*c_e;    
    fp(Np+1) = -2.0*kappaL/hp*(phi2(i+iAdd) - phi2_p((i-1)*Np + Np)); % current flux
    fc(Np+1) = -2.0*DMes/hp*(c(i+iAdd) - cp((i-1)*Np + Np)); % mass flux
    
    % particle (compute balances)
    for j=1:Np
        ddphi_pdt((i-1)*Np + j) = -1/aCMes(j)*(fp(j)*Sr(j) - fp(j+1)*Sr(j+1))/V(j);
        
        dcpdt((i-1)*Np + j) = 1.0/epsMes*(fc(j)*Sr(j) - fc(j+1)*Sr(j+1))/V(j) + ...
            aCMes(j)/2/epsMes/F*ddphi_pdt((i-1)*Np + j);
    end
    
end

% % supress temporarily the equations for potential
% dphidt(:) = 0.0; 
% dphi_pdt(:) = 0.0;

% % supress temporarily the equations for concentration
% dcdt(:) = 0.0; 
% dcpdt(:) = 0.0;

dydt = [dphi2dt; ddphidt; dcdt; ddphi_pdt; dcpdt];

end

























