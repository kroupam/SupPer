function EIS = EIS(par)
% This file is part of the program SupPer for the modelling of
% dynamic behavior and performance of Supercapacitors
% (Electrochemical Double Layer Capacitors)
% 
% Author: Martin Kroupa
% Oct 2015 - Mar 2016
% 
% Developed jointly at Imperial College London and University of Chemistry
% and Technology Prague

disp('Entering the Electrochemical Impedance Spectroscopy calculation')

%% Unpack parameters
v2struct(par.spec); % unpack specification parameters
v2struct(par.gen); % unpack general parameters
v2struct(par.num); % unpack numerical parameters
v2struct(par.EIS); % unpack EIS parameters

switch dim
    case '1D1D'
        v2struct(par.mes); % unpack parameters of meso-pores
end

% initialize fields
Z = NaN*zeros(N_f,1);
Cs = NaN*zeros(N_f,1);
CsIm = NaN*zeros(N_f,1);
lambda = NaN*zeros(N_f,1);

for k=1:N_f
    
    fprintf('\n'); % new line in output
    disp(['Computing frequency f = ' num2str(fVec(k)) ' 1/s'])
    
    % Intitial Conditions -----------------------------------------------------
    phiInit = (UFix + UAmp)/2; % initial condition (potential in macro-pores)
    p20 = ones(2*N+Ns,1);
    p20 = p20*phiInit; % consistent for U0 = phiInit
    
    dp0 = ones(2*N,1); % initial condition for solid phase potential
    dp0(1:N) = -dp0(1:N)*phiInit; % left electrode
    dp0(N+1:2*N) = dp0(N+1:2*N)*phiInit; % right electrode
    
    % initial condition (concentration in macro-pores)
    c0 = ones(2*N+Ns,1)*cB;
    
    switch dim
        case '1D'
            y0 = [p20; dp0; c0]; % pack initial conditions
        case '1D1D'
            % initial condition (potential in meso-pores)
            dpp0 = ones(2*N*Np,1);
            dpp0(1:N*Np) = -dpp0(1:N*Np)*phiInit; % left electrode
            dpp0(N*Np+1:2*N*Np) = dpp0(N*Np+1:2*N*Np)*phiInit; % right electrode
            
            % initial condition (concentration in meso-pores)
            cp0 = ones(2*N*Np,1)*cB;
            
            y0 = [p20; dp0; c0; dpp0; cp0]; % pack initial conditions
    end
    
    tStep = 1.0/fVec(k)/resT; % time-step [s]
    tBeg = 0.0; % start time of simulation [s]
    tEnd = nT/fVec(k); % end time of simulation [s]
    
    options = odeset( ...
        'MaxStep',tStep, ...
        'RelTol', rTol, ...
        'AbsTol', aTol, ...
        'MassSingular', 'yes', ...
        'MStateDependence', 'none', ...
        'Mass', M ...
        );
    
    % wrap parameters
    parIn = fVec(k);
    
    % Integrate in time -------------------------------------------------------
    [tt,yy] = ode15s(@(t,y) rhsEIS(t,y,par,parIn), tBeg:tStep:tEnd, y0, options);
    
    phi2 = yy(:,1:2*N+Ns); % potential in electrolyte (phi2, algebraic variable)
    dphi = yy(:,2*N+Ns+1:4*N+Ns); % potential difference (phi1-phi2)
    c = yy(:,4*N+Ns+1:6*N+2*Ns); % concentration
    
    % extract phi1 (potential in the solid phase)
    phi1 = NaN*phi2;
    phi1(:,1:N) = dphi(:,1:N) + phi2(:,1:N);
    phi1(:,N+Ns+1:2*N+Ns) = dphi(:,N+1:2*N) + phi2(:,N+Ns+1:2*N+Ns);

    UU = UFix + UAmp*exp(1i*2*pi*fVec(k)*tt); % compute voltage
    jj = 2.0*sigmaS*phi1(:,1)/h; % current density at the left current collector
   
    idx = (nT-1)*resT; % the beginning of the last period
    
    I_guess = 10*Cmax*UAmp*fVec(k); % initial guess of the current amplitude
    
%     disp(['I_guess = ' num2str(I_guess)])
    
    % fit the phase angle from the dependence of current on time
    [x,~,exitflag] = lsqnonlin(@(cons) impFit(cons,tt(idx:end),jj(idx:end)*S,fVec(k)), ...
        [I_guess 0.5]);
    
    if exitflag <=0
        error('Failed to fit impedance')
    end
    
    j0 = x(1); % extract current amplitude
    phi_s = x(2); % extract phase angle
    
    disp(['j0 = ' num2str(j0)])
    disp(['phi_s = ' num2str(phi_s)])
    
    Z(k) = UAmp/j0*exp(-1i*phi_s); % compute complex impedance
    
    Za = sqrt(real(Z(k))^2 + imag(Z(k))^2); % compute impedance modulus
    
    Cs(k) = -imag(Z(k))/(2*pi*fVec(k))/Za^2; % compute capacitance

    CsIm(k) = real(Z(k))/(2*pi*fVec(k))/Za^2; % imaginary part of capacitance
    
    phiLo = min(real(dphi(idx:end,1:2*N)),[],1); % minimum potential during simulation
    phiHi = max(real(dphi(idx:end,1:2*N)),[],1); % maximum potential during simulation
    
    sumL = 1/(N-1)*trapz(phiHi(1:N)-phiLo(1:N)); % summation
    sumR = 1/(N-1)*trapz(phiHi(N+1:2*N)-phiLo(N+1:2*N)); % summation
    
    lambda(k) = (sumL+sumR)/UAmp/2; % penetration depth

    disp(['C = ' num2str(Cs(k)) ' F'])
    disp(['lambda = ' num2str(lambda(k))])
    
    cMin = min(c,[],2); % get the minumum concentration in the system
    cMax = max(c,[],2); % get the maximum concentration in the system
    
    switch dim
        case '1D1D'
%             dphi_p = yy(:,6*N+2*Ns + 1:6*N+2*Ns + 2*N*Np);  % potential difference (phi1-phi2p) in meso-pores
            cp = yy(:,6*N+2*Ns + 2*N*Np + 1:6*N+2*Ns + 4*N*Np);  % concentration of electrolyte in meso-pores
            
            cpMin = min(cp,[],2); % get the minumum concentration in particles
            cpMax = max(cp,[],2); % get the maximum concentration in particles
            
            cMin = min(cMin,cpMin);  % get minimum of both
            cMax = max(cMax,cpMax);  % get maximum of both
            
            phi2_p = NaN*dphi_p;
            for i=1:N
                for j=1:length(tt)
                    phi2_p(j,(i-1)*Np+1:(i-1)*Np+Np) = phi1(j,i) - ...
                        dphi_p(j,(i-1)*Np+1:(i-1)*Np+Np);
                end
            end
            for i=N+1:2*N
                for j=1:length(tt)
                    phi2_p(j,(i-1)*Np+1:(i-1)*Np+Np) = phi1(j,Ns+i) - ...
                        dphi_p(j,(i-1)*Np+1:(i-1)*Np+Np);
                end
            end
    end
    
    % store this cycle in the overall vectors
    EIS.f{k}.tt = tt;
    EIS.f{k}.phi2 = phi2;
    EIS.f{k}.phi1 = phi1;
    EIS.f{k}.c = c;
    EIS.f{k}.jj = jj;
    EIS.f{k}.UU = UU;
    EIS.f{k}.Cs = Cs(k);
    EIS.f{k}.CsIm = CsIm(k);
    EIS.f{k}.Z = Z(k);
    EIS.f{k}.lambda = lambda(k);
    EIS.f{k}.cMin = cMin;
    EIS.f{k}.cMax = cMax;
    
    switch dim
        case '1D1D'
            EIS.f{k}.phi2_p = phi2_p;
            EIS.f{k}.cp = cp;
    end
    
end

% compute time constant
td = NaN; % reset
Cs_max = max(Cs);
for k=1:N_f
    if Cs(k) <= Cs_max/2
        if k==1 % something wrong probably happened
            td = NaN; % time constant
            break
        end
            
        % interpolate knee frequency
        fK = fVec(k-1) + (fVec(k) - fVec(k-1))* ...
            (Cs_max/2 - Cs(k-1))/(Cs(k) - Cs(k-1)); % lever rule
        td = 1.0/fK; % time constant
        break
    end
end

EIS.td = td;






