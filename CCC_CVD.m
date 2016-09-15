clear all; close all; clc %#ok<CLSCR>
% This file is part of the program SupPer for the modelling of
% dynamic behavior and performance of Supercapacitors
% (Electrochemical Double Layer Capacitors)
% 
% Author: Martin Kroupa
% Oct 2015 - Mar 2016
% 
% Developed jointly at Imperial College London and University of Chemistry
% and Technology Prague

% Standalone program to compute Constant Current Charging followed by
% Constant Voltage Discharging for comparison with literature data

%% Initialize (load default parameters)
disp('Starting: Constant Current Charging followed by Constant Voltage Discharging');

aa_Ver = 4.5e7;% specific surface area [m^2/m^3] (to fit Verbrugge)
% aa_Ver = 4.8e7;% specific surface area [m^2/m^3] (to fit Verbrugge)

par.spec = parSpec; % load specifications about what should be computed
par.mes = []; % assign dummy field, only to satisfy format
par.num = []; % assign dummy field, only to satisfy format
par.gen = []; % assign dummy field, only to satisfy format
par.gen.aaTot = aa_Ver;
par = parameters(par);
par.CCC = parCCC; % load parameters for the operation mode
v2struct(par.spec);
v2struct(par.gen);
v2struct(par.num);
v2struct(par.CCC);

tt = [];
jj = [];
UU = [];
phi2 = [];
phi1 = [];
c = [];
cMin = [];
cMax = [];

%% constant current charging

js = I_fix/S; % current density [A/m^2]

% Intitial Conditions -----------------------------------------------------
% initial condition (potential in macro-pores)

kappaL = F^2/R/T*2*DS*cB;
d_phi_fix = -js*(Ls/2)./kappaL; % fixed potetnial difference accros separator [V]
% d_phi_fix = 0.0;
phi_Init_L = UMin/2; % + d_phi_fix; % Initial condition in the left electrode
p20 = ones(2*N+Ns,1);
p20 = p20*phi_Init_L;

dp0 = ones(2*N,1); % initial condition for solid phase potential
dp0(1:N) = -dp0(1:N)*(UMin/2 + d_phi_fix);
dp0(N+1:2*N) = dp0(N+1:2*N)*(UMin/2 + d_phi_fix);

% initial condition (concentration in macro-pores)
c0 = ones(2*N+Ns,1)*cB;
y0 = [p20; dp0; c0];

tBeg = 0.0; % start time of simulation [s]
tEnd = 17.9; % end time of simulation [s]
% tEnd = 23.2; % end time of simulation [s]
tStep = tEnd/2.0/resT; % time-step [s]

% wrap parameters
parIn = js;

options = odeset( ...
    'MaxStep',tStep, ...
    'RelTol', rTol, ...
    'AbsTol', aTol, ...
    'MassSingular', 'yes', ...
    'MStateDependence', 'none', ...
    'Mass', M ...
    );

% Integrate in time -------------------------------------------------------
[tt_t,yy] = ode15s(@(t,y) rhsCCC(t,y,par,parIn), tBeg:tStep:tEnd, y0, options);

phi2_t = yy(:,1:2*N+Ns); % potential in electrolyte (phi2, algebraic variable)
phiDelta_t = yy(:,2*N+Ns+1:4*N+Ns); % potential difference (phi1-phi2)
c_t = yy(:,4*N+Ns+1:6*N+2*Ns); % concentration

% extract phi1
phi1_t = NaN*phi2_t;
phi1_t(:,1:N) = phiDelta_t(:,1:N) + phi2_t(:,1:N);
phi1_t(:,N+Ns+1:2*N+Ns) = phiDelta_t(:,N+1:2*N) + phi2_t(:,N+Ns+1:2*N+Ns);

U_t = js*h/2.0/sigmaS + phi1_t(:,2*N+Ns);
j_t = 2.0*sigmaS*phi1_t(:,1)/h; % current density at the left current collector

% kappaL = F^2/R/T*2*DS*cB;
% d_phi_fix = -js*hs/2.0/kappaL; % fixed potetnial difference accros separator [V]
% U_t = 2*(phi_t(:,N+Ns) - d_phi_fix);
% % j_t = kappaL.*(U_t(:,1)/2 - phi_t(:,N))/(h/2+Ls/2);
% j_t = 2.0*kappaL.*(U_t(:,1)/2 - phi_t(:,N+Ns))/hs;

cMin_t = min(c_t,[],2); % get the minumum concentration in the system
cMax_t = max(c_t,[],2); % get the maximum concentration in the system

tt = [tt; tt_t(:)];
jj = [jj; j_t(:)];
UU = [UU; U_t(:)];
phi2 = [phi2; phi2_t(:,:)];
phi1 = [phi1; phi1_t(:,:)];
c = [c; c_t(:,:)];
cMin = [cMin; cMin_t(:)];
cMax = [cMax; cMax_t(:)];

%% constant voltage discharging

% initial condition
p20(:) = phi2_t(end,:)';
dp0(:) = phiDelta_t(end,:)';
c0(:) = c_t(end,:)';
y0 = [p20; dp0; c0];

tBeg = tEnd; % start time of simulation [s]
tEnd = 30.0; % end time of simulation [s]

% wrap parameters
parIn = UDis;

% Integrate in time -------------------------------------------------------
[tt_t,yy] = ode15s(@(t,y) rhsCVD(t,y,par,parIn), tBeg:tStep:tEnd, y0, options);

phi2_t = yy(:,1:2*N+Ns); % potential in electrolyte (phi2, algebraic variable)
phiDelta_t = yy(:,2*N+Ns+1:4*N+Ns); % potential difference (phi1-phi2)
c_t = yy(:,4*N+Ns+1:6*N+2*Ns); % concentration

% extract phi1
phi1_t = NaN*phi2_t;
phi1_t(:,1:N) = phiDelta_t(:,1:N) + phi2_t(:,1:N);
phi1_t(:,N+Ns+1:2*N+Ns) = phiDelta_t(:,N+1:2*N) + phi2_t(:,N+Ns+1:2*N+Ns);

U_t = UDis*ones(length(tt_t),1);
j_t = 2.0*sigmaS*phi1_t(:,1)/h; % current density at the left current collector

% kappaL = F^2/R/T*2*DS*cB;
% % j_t = kappaL.*(U_t(:,1)/2 - phi_t(:,N))/(h/2+Ls/2);
% j_t = 2.0*kappaL.*(U_t(:,1)/2 - phi_t(:,N+Ns))/hs;

cMin_t = min(c_t,[],2); % get the minumum concentration in the system
cMax_t = max(c_t,[],2); % get the maximum concentration in the system

kappaL = F^2/R/T*2*D*cB;
kappaLS = F^2/R/T*2*DS*cB;
RR_Liq = 2*Le/kappaL/S;
RR_Sol = 2*Le/sigmaS/S;
RR_E = (RR_Liq + RR_Sol)/2;
RR_S = Ls/kappaLS/S;
RR = RR_E + RR_S;

% RR = (Le)/kappaL/S;

tau = RR*Cmax;

I_Charg = 100; % [A]

V_RC1 = 1.6 + I_Charg/Cmax*tt;
V_RC1(2:end) = V_RC1(2:end) + (RR_S+RR_Sol)*I_Charg;

ttRC = tt_t - min(tt_t);
V_RC = (V_RC1(end) - U_t(end))*exp(-ttRC/tau);
I_RC = V_RC/RR;

tt = [tt; tt_t(:)];
jj = [jj; j_t(:)];
UU = [UU; U_t(:)];
phi2 = [phi2; phi2_t(:,:)];
phi1 = [phi1; phi1_t(:,:)];
c = [c; c_t(:,:)];
cMin = [cMin; cMin_t(:)];
cMax = [cMax; cMax_t(:)];

res.tt = tt;
res.UU = UU;
res.jj = jj;
res.phi2 = phi2;
res.phi1 = phi1;
res.c = c;
res.cMin = cMin;
res.cMax = cMax;

%% Plot spatial distributions

k = 1; % index of frequency / scan rate to plot
% tBeg = 17.9; % initial time
% tEnd = 30.0; % end time
tBeg = 0.0; % initial time
tEnd = 17.8; % end time
par.spec.method = 'CCC';

% plot results
figurePotConXTime(res, par, k, tBeg, tEnd)

% print('-f1','-dpng','figures/final/pot2XTimeBeg')
% print('-f2','-dpng','figures/final/pot1XTimeBeg')
% print('-f3','-dpng','figures/final/conXTimeBeg')

% print('-f1','-dpng','figures/final/pot2XTimeEnd')
% print('-f2','-dpng','figures/final/pot1XTimeEnd')
% print('-f3','-dpng','figures/final/conXTimeEnd')


%% Plot current and voltage characteristics

figureCurVolTimeVerbrugge(res, par,tt_t, I_RC, V_RC1)

% print('-f1','-dpng','figures/final/Verb_U')
% print('-f2','-dpng','figures/final/Verb_I')








