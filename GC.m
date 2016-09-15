function GC = GC(par)
% This file is part of the program SupPer for the modelling of
% dynamic behavior and performance of Supercapacitors
% (Electrochemical Double Layer Capacitors)
% 
% Author: Martin Kroupa
% Oct 2015 - Mar 2016
% 
% Developed jointly at Imperial College London and University of Chemistry
% and Technology Prague

disp('Entering the Galvanostatic Cycling calculation')

%% Unpack parameters
v2struct(par.spec); % unpack specification parameters
v2struct(par.gen); % unpack general parameters
v2struct(par.num); % unpack numerical parameters
v2struct(par.GC); % unpack GC parameters

switch dim
    case '1D1D'
        v2struct(par.mes); % unpack parameters of meso-pores
end

Cs = NaN*zeros(N_I,1); % initialize capacitance
for k=1:N_I
    skip = 0; % flag determining whether the integration is to be skipped
    
    fprintf('\n'); % new line in output
    disp(['Computing current I = ' num2str(IVec(k)) ' A'])

    % Intitial Conditions -----------------------------------------------------
    kappaL = F^2/R/T*2*DS*cB;
    js = IVec(k)/S; % current density [A/m^2]
    U0 = -js*Ls./kappaL; % Initial Voltage increase due to fixed potetnial difference accros separator [V]
    
    phiInit = UMin/2; % initial condition (potential in macro-pores)
    p20 = ones(2*N+Ns,1);
    p20 = p20*phiInit;
    
    dp0 = ones(2*N,1); % initial condition for solid phase potential
    dp0(1:N) = -dp0(1:N)*UMin/2; % left electrode
    dp0(N+1:2*N) = dp0(N+1:2*N)*UMin/2; % right electrode
        
    % initial condition (concentration in macro-pores)
    c0 = ones(2*N+Ns,1)*cB;
    
    if abs(U0) >= abs(UMax - UMin) % the Ohmic drop is larger than the voltage window
        skip = 1;
    end
    
    switch dim
        case '1D1D'
            % initial condition (potential in meso-pores)
            dpp0 = ones(2*N*Np,1);
            dpp0(1:N*Np) = -dpp0(1:N*Np)*phiInit; % left electrode
            dpp0(N*Np+1:2*N*Np) = dpp0(N*Np+1:2*N*Np)*phiInit; % right electrode
            
            % initial condition (concentration in meso-pores)
            cp0 = ones(2*N*Np,1)*cB;
            
            phi2_p = [];
            cp = [];
    end
    
    tSim = 2.0*Cmax*abs(UMax - UMin)/IVec(k); % estimated time of simulation (safe guess)
    
    tStep = tSim/2.0/resT; % time-step [s]
    tEnd = 0.0; % start time of simulation [s]
    
    tt = [];
    phi2 = [];
    phi1 = [];
    c = [];
    jj = [];
    UU = [];
    cMin = [];
    cMax = [];
    Cs_n = NaN*zeros(NC,1);
    lambda = NaN*zeros(NC,1);
    
    flag = 0;
    for n = 1:NC
        if skip == 1 % the integration would be meaningless => skip it
            disp('Ohmic drop is larger than voltage window, C = NaN')
            break
        end
        disp(['cycle number n = ' int2str(n)])
        
        tBeg0 = tEnd; % start time of the cycle
        
        for p=1:2 % charging and discharging for each cycle
            
            switch p
                case 1
                    js = IVec(k)/S; % current density [A/m^2]
                case 2
                    js = -IVec(k)/S; % current density [A/m^2]
            end
            
            % wrap parameters
            parIn = js;
            
            options = odeset( ...
                'MaxStep',tStep, ...
                'RelTol', rTol, ...
                'AbsTol', aTol, ...
                'Events', @(t,y) eventGC(t,y,par,parIn), ...
                'MassSingular', 'yes', ...
                'MStateDependence', 'none', ...
                'Mass', M ...
                );
            
            switch dim
                case '1D'
                    y0 = [p20; dp0; c0]; % pack initial conditions
                case '1D1D'
                    y0 = [p20; dp0; c0; dpp0; cp0]; % pack initial conditions
            end
            
            tBeg = tEnd; % start time of simulation [s]
            tEnd = tBeg + tSim; % end time of simulation [s]
            
            % Integrate in time -------------------------------------------------------
            [tt_t,yy,tE_out] = ode15s(@(t,y) rhsGC(t,y,par,parIn), tBeg:tStep:tEnd, y0, options);
            
            if isempty(tE_out) || length(tE_out) > 1 || tE_out ~= tt_t(end) % the ohmic dop is larger than U => break
                flag = 1; % set flag to break from the upper loop
                tEnd = NaN;
                disp('Integration failed, C = NaN')
                break
            else
                tEnd = tE_out;
            end
            
            phi2_t = yy(:,1:2*N+Ns); % potential in electrolyte (phi2, algebraic variable)
            dphi_t = yy(:,2*N+Ns+1:4*N+Ns); % potential difference (phi1-phi2)
            c_t = yy(:,4*N+Ns+1:6*N+2*Ns); % concentration
            
            % extract phi1 (potential in the solid phase)
            phi1_t = NaN*phi2_t;
            phi1_t(:,1:N) = dphi_t(:,1:N) + phi2_t(:,1:N);
            phi1_t(:,N+Ns+1:2*N+Ns) = dphi_t(:,N+1:2*N) + phi2_t(:,N+Ns+1:2*N+Ns);
            
            U_t = js*h/2.0/sigmaS + phi1_t(:,2*N+Ns); % compute voltage
            j_t = 2.0*sigmaS*phi1_t(:,1)/h; % current density at the left current collector
            
            phiLo = min(dphi_t(:,1:2*N),[],1); % minimum potential during simulation
            phiHi = max(dphi_t(:,1:2*N),[],1); % maximum potential during simulation
            
            sumL = 1/(N-1)*trapz(phiHi(1:N)-phiLo(1:N)); % summation
            sumR = 1/(N-1)*trapz(phiHi(N+1:2*N)-phiLo(N+1:2*N)); % summation
            
            lambda(n) = (sumL+sumR)/(UMax-UMin); % penetration depth
            
            % assign intitial conditions for next cycle
            p20 = phi2_t(end,:)';
            dp0 = dphi_t(end,:)';
            c0 = c_t(end,:)';
            
            cMin_t = min(c_t,[],2); % get the minumum concentration in the system
            cMax_t = max(c_t,[],2); % get the maximum concentration in the system
            
            switch dim
                case '1D1D'
                    dphi_p_t = yy(:,6*N+2*Ns + 1:6*N+2*Ns + 2*N*Np); % potential difference (phi1-phi2p) in meso-pores
                    cp_t = yy(:,6*N+2*Ns + 2*N*Np + 1:6*N+2*Ns + 4*N*Np); % concentration of electrolyte in meso-pores
                    
                    % assign intitial conditions for next cycle
                    dpp0 = dphi_p_t(end,:)';
                    cp0 = cp_t(end,:)';
                    
                    cpMin_t = min(cp_t,[],2); % get the minumum concentration in particles
                    cpMax_t = max(cp_t,[],2); % get the maximum concentration in particles
                    
                    cMin_t = min(cMin_t,cpMin_t);  % get minimum of both
                    cMax_t = max(cMax_t,cpMax_t);  % get maximum of both
                    
                    phi2_p_t = NaN*dphi_p_t;
                    for i=1:N
                        for j=1:length(tt_t)
                            phi2_p_t(j,(i-1)*Np+1:(i-1)*Np+Np) = phi1_t(j,i) - ...
                                dphi_p_t(j,(i-1)*Np+1:(i-1)*Np+Np);
                        end
                    end
                    for i=N+1:2*N
                        for j=1:length(tt_t)
                            phi2_p_t(j,(i-1)*Np+1:(i-1)*Np+Np) = phi1_t(j,Ns+i) - ...
                                dphi_p_t(j,(i-1)*Np+1:(i-1)*Np+Np);
                        end
                    end
                    
                    phi2_p = [phi2_p; phi2_p_t(:,:)];
                    cp = [cp; cp_t(:,:)];
            end
            
            % store this cycle in the overall vectors
            tt = [tt; tt_t(:)];
            phi2 = [phi2; phi2_t(:,:)];
            phi1 = [phi1; phi1_t(:,:)];
            c = [c; c_t(:,:)];
            jj = [jj; j_t(:)];
            UU = [UU; U_t(:)];
            cMin = [cMin; cMin_t(:)];
            cMax = [cMax; cMax_t(:)];
            
        end
        
        if flag == 1 % the integration failed in the inner loop
            break
        end

        % compute capacity from one cycle
        Cs_n(n) = IVec(k)*(tEnd-tBeg0)/2/(UMax - UMin);
        disp(['C = ' num2str(Cs_n(n)) ' F'])
        disp(['lambda = ' num2str(lambda(n))])
        
        % conditional break if steady-state capacitance is reached
        if n~=1
            if abs(Cs_n(n) - Cs_n(n-1)) < thresh*Cmax
                % shorten vector
                Cs_n = Cs_n(1:n);
                lambda = lambda(1:n);
                break
            end
        end
        
    end
    
    Cs(k) = Cs_n(end);

    freq = IVec(k)/Cmax/(2.0*(UMax - UMin)); % compute frequency
    
    % save fields to exported variables (structure res)
    GC.I{k}.tt = tt;
    GC.I{k}.phi2 = phi2;
    GC.I{k}.phi1 = phi1;
    GC.I{k}.c = c;
    GC.I{k}.jj = jj;
    GC.I{k}.UU = UU;
    GC.I{k}.Cs = Cs(k);
    GC.I{k}.f = freq; % one hypothetical cycle
    GC.I{k}.lambda = lambda(end);
    GC.I{k}.cMin = cMin;
    GC.I{k}.cMax = cMax;
    
    switch dim
        case '1D1D'
            GC.I{k}.phi2_p = phi2_p;
            GC.I{k}.cp = cp;
    end
    
end

% compute time constant
td = NaN; % reset
Cs_max = max(Cs);
for k=1:N_I
    if Cs(k) < Cs_max/2
        % interpolate knee frequency
        fK = GC.I{k-1}.f + (GC.I{k}.f - GC.I{k-1}.f)* ...
            (Cs_max/2 - Cs(k-1))/(Cs(k) - Cs(k-1)); % lever rule
        td = 1.0/fK; % time constant
        break
    end
end

GC.td = td;




















