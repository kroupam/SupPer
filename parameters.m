%% Load all default parameters
% This file is part of the program SupPer for the modelling of
% dynamic behavior and performance of Supercapacitors
% (Electrochemical Double Layer Capacitors)
% 
% Author: Martin Kroupa
% Oct 2015 - Mar 2016
% 
% Developed jointly at Imperial College London and University of Chemistry
% and Technology Prague

function par = parameters(par)

v2struct(par.spec) % unpack specification parameters

par.fig = parFigure; % load parameters for figures
par.gen = parGeneral(par); % load general physical parameters and constants
switch dim
    case '1D1D'
        par.mes = parMeso(par); % load parameters of mesopores
end
par.num = parNumerical(par); % load numerical parameters

end

%% General physical parameters and constants
function gen = parGeneral(par)

% Physical constants 
e = 1.602e-19; % elementary charge [C]
NA = 6.022e23; % Avogadro number [mol^-1]
R = 8.314; % universal gas constant [J/K/mol]
perm0 = 8.854e-12; % vacuum perimttivity [F/m]
F = e*NA; % Faraday constant [C/mol]
kB = R/NA; % Boltzmann constant [m^2*kg/s^2/K]

% Geometrical parameters
Le = 50.0e-6; % width of electrode [m]
S = 2.747; % area of the electrode [m^2]
Ls = 25.0e-6; % width of separator [m]

% Parameters of electrode
epsTot = 0.67; % total porosity of electrode
gamTot = 2.3; % tortuosity of macro-pores
aaTot = 3.89e7;% specific surface area [m^2/m^3]
sigmaS = 0.0521; % solid-phase conductivity [S/m]

% Parameters of electrolyte
D0 = 9.6e-12; % Diffusion coefficient in macro-pores [m^2/s] (uncorrected)
cB = 0.93*1e3; % bulk concentration [mol/m^3]
T = 298; % temperature [K]
D = D0*epsTot/gamTot; % Corrected diffusion coefficient in macro-pores [m^2/s]
permR = 36.6; % relative permittivity (Acetonitrile (Sigma Aldrich))
perm = perm0*permR; % total permittivity
lambdaS = 0.30e-9; % thickness of Stern layer [m] (half of bare ion size)
a = 1.2e-9; % effective solvated ion diameter [m] (for Bikerman model)

Cmax = aaTot*perm/lambdaS*S*Le/2; % maximum capacity

% Parameters of separator
epsS = 0.6; % porosity of separator
gamS = 1.29; % tortuosity of separator
DS = D0*epsS/gamS; % Corrected diffusion coefficient in separator [m^2/s]

% Pack general parameters
listGen = {'fieldNames' ...
    ,'e' ...
    ,'R' ...
    ,'F' ...
    ,'NA' ...
    ,'R' ...
    ,'Le' ...
    ,'S' ...
    ,'Ls' ...
    ,'epsTot' ...
    ,'gamTot' ...
    ,'aaTot' ...
    ,'Cmax' ...
    ,'cB' ...
    ,'T' ...
    ,'D0' ...
    ,'D' ...
    ,'perm' ...
    ,'lambdaS' ...
    ,'a' ...
    ,'epsS' ...
    ,'DS' ...
    ,'sigmaS' ...
    };

gen = v2struct(listGen);
end

%% Parameters of meso pores
function mes = parMeso(par)

% unpack general parameters
epsTot = par.gen.epsTot;
D0 = par.gen.D0;
aaTot = par.gen.aaTot;

% Parameters of electrode
epsMac = 0.36; % macro-porosity (must be smaller than eps)
if epsTot<epsMac
    error('epsTot must be larger than epsMac')
end
epsMes = epsTot - epsMac; % meso-porosity
epsMes = epsMes/(1.0-epsMac); % Transform to the framework of particles

gamMac = 1.4; % tortuosity of macro-pores
gamMes = 2.6; % tortuosity of meso-pores

% Parameters of electrolyte
DMac = D0*epsMac/gamMac; % Corrected diffusion coefficient in macro-pores [m^2/s]
DMes = D0*epsMes/gamMes; % Corrected diffusion coefficient in meso-pores [m^2/s]

% Parameters of meso-pores
Rp = 500.0e-9; % radius of carbon particles [m]

aaMac = 3.0*(1.0-epsMac)/Rp; % specific surface area of macro-pores [m^2/m^3]

if aaTot < aaMac
    error('aaTot must be larger than aaMac')
end

aaMes = aaTot - aaMac; % specific surface area of meso-pores [m^2/m^3]

aaMes = aaMes/(1.0-epsMac); % Transform to the framework of particles

% Pack meso parameters
listMes = {'fieldNames' ...
    ,'epsMac' ...
    ,'epsMes' ...
    ,'DMac' ...
    ,'DMes' ...
    ,'Rp' ...
    ,'aaMes' ...
    ,'aaMac' ...
    };

mes = v2struct(listMes);
end

%% Numerical Parameters
function num = parNumerical(par)

% Extract necessary parameters
Le = par.gen.Le;
Ls = par.gen.Ls;
dim = par.spec.dim;

% time-resolution (used to compute time-step)
resT = 500; % resolution of the time output (default value)

% Discretization parameters
N = 30; % number of control volumes in electrode in x-direction
Ns = 15; % number of control volumes in separator

h = Le/N; % discretization step in x-direction
hs = Ls/Ns; % discretization step in separator

% Parameters of integration
rTol = 1e-6; % relative tolerance

switch dim
    case '1D'
        aTol = zeros(6*N+2*Ns,1);
        aTol(1:2*N+Ns,1) = 1.0e-7; % potential in electrolyte (phi2: algebraic)
        aTol(2*N+Ns+1:4*N+Ns,1) = 1.0e-7; % potential difference (phi1-phi2)
        aTol(4*N+Ns+1:6*N+2*Ns,1) = 1.0e-3; % concentration
        
        M = eye(6*N+2*Ns); % create mass matrix of the DAE system
        
        listNum = {'fieldNames'}; % initialize
    case '1D1D'
        Rp = par.mes.Rp;
        % Discretization parameters for meso-pores
        Np = 20; % number of control volumes in y-direction
        hp = Rp/Np; % discretization step in the r-direction
        
        r = zeros(Np+1,1);
        Sr = zeros(Np+1,1);
        V = zeros(Np,1);
        
        for i=1:Np+1
            r(i) = (i-1)*hp; % radial position of the boundary
            Sr(i) = 4.0*pi*r(i)^2; % area between control volumes
        end
        
        for i=1:Np
            V(i) = 4.0/3.0*pi*(r(i+1)^3-r(i)^3); % volume of control volume
        end
        
        aTol = zeros(6*N+2*Ns+4*N*Np,1);
        aTol(1:2*N+Ns,1) = 1.0e-7; % potential in electrolyte in macro-pores (phi2: algebraic)
        aTol(2*N+Ns+1:4*N+Ns,1) = 1.0e-7; % potential difference (phi1-phi2)
        aTol(4*N+Ns+1:6*N+2*Ns,1) = 1.0e-3; % concentration in macro-pores
        aTol(6*N+2*Ns + 1:6*N+2*Ns + 2*N*Np) = 1.0e-7; % potential difference in meso-pores (phi1-phi2p)
        aTol(6*N+2*Ns + 2*N*Np + 1:6*N+2*Ns + 4*N*Np) = 1.0e-3; % concentration in meso-pores
        
        M = eye(6*N+2*Ns+4*N*Np); % create mass matrix of the DAE system
        
        listNum = {'fieldNames' ...
            ,'Np' ...
            ,'hp' ...
            ,'Sr' ...
            ,'V' ...
            };
end

for i=1:2*N+Ns
    M(i,i) = 0; % nulify positions with algebraic variables (potential electrolyte)
end

thresh = 1.0e-4; % maximum relative change of capacity to be considered as a steady state

len = length(listNum);
listNum{len+1} = 'N';
listNum{len+2} = 'h';
listNum{len+3} = 'rTol';
listNum{len+4} = 'aTol';
listNum{len+5} = 'thresh';
listNum{len+6} = 'resT';
listNum{len+7} = 'Ns';
listNum{len+8} = 'hs';
listNum{len+9} = 'M';

% Pack numerical parameters
num = v2struct(listNum);
end













































