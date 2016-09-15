function CV = CV(par)
% This file is part of the program SupPer for the modelling of
% dynamic behavior and performance of Supercapacitors
% (Electrochemical Double Layer Capacitors)
% 
% Author: Martin Kroupa
% Oct 2015 - Mar 2016
% 
% Developed jointly at Imperial College London and University of Chemistry
% and Technology Prague

disp('Entering the Cyclic Voltammetry calculation')

%% Unpack parameters
v2struct(par.spec); % unpack specification parameters
v2struct(par.gen); % unpack general parameters
v2struct(par.num); % unpack numerical parameters
v2struct(par.CV); % unpack CV parameters

switch dim
    case '1D1D'
        v2struct(par.mes); % unpack parameters of meso-pores
end

n_t = 2*resT + 1; % number of time outputs

Cs = NaN*zeros(N_nu,1); % initialize capacitance
for k=1:N_nu
    
    fprintf('\n'); % new line in output
    disp(['Computing scan rate nu = ' num2str(nuVec(k)) ' V/s'])
    
    t0 = (UMax - UMin)/nuVec(k); % half cycle period

    % Intitial Conditions -----------------------------------------------------
    phiInit = UMin/2; % initial condition (potential in macro-pores)
    p20 = ones(2*N+Ns,1);
    p20 = p20*phiInit; % consistent for U0 = 2*phiInit
    
    dp0 = ones(2*N,1); % initial condition for solid phase potential
    dp0(1:N) = -dp0(1:N)*UMin/2; % left electrode
    dp0(N+1:2*N) = dp0(N+1:2*N)*UMin/2; % right electrode

    % initial condition (concentration in macro-pores)
    c0 = ones(2*N+Ns,1)*cB;
    
    switch dim
        case '1D1D'
            % initial condition (potential in meso-pores)
            dpp0 = ones(2*N*Np,1);
            dpp0(1:N*Np) = -dpp0(1:N*Np)*phiInit; % left electrode
            dpp0(N*Np+1:2*N*Np) = dpp0(N*Np+1:2*N*Np)*phiInit; % right electrode
            
            % initial condition (concentration in meso-pores)
            cp0 = ones(2*N*Np,1)*cB;
            
%             phi2_p = zeros(NC*n_t,2*N*Np);
%             cp = zeros(NC*n_t,2*N*Np);
    end
    
    tStep = (UMax - UMin)/nuVec(k)/resT; % time-step [s]
    tEnd = 0.0; % start time of simulation [s]
    
    options = odeset( ...
        'MaxStep',tStep ...
        ,'RelTol', rTol ...
        ,'AbsTol', aTol ...
        ,'MassSingular', 'yes' ...
        ,'MStateDependence', 'none' ...
        ,'Mass', M ...
        );
    
    tt = zeros(NC*n_t,1);
    phi2 = zeros(NC*n_t,2*N+Ns);
    phi1 = zeros(NC*n_t,2*N+Ns);
    c = zeros(NC*n_t,2*N+Ns);
    jj = zeros(NC*n_t,1);
    UU = zeros(NC*n_t,1);
    cMin = zeros(NC*n_t,1);
    cMax = zeros(NC*n_t,1);
    Cs_n = zeros(NC,1);
    lambda = zeros(NC,1);
    
    for n = 1:NC        
        disp(['cycle number n = ' int2str(n)])
        
        % wrap parameters
        parIn = nuVec(k);
        
        switch dim
            case '1D'
                y0 = [p20; dp0; c0]; % pack initial conditions
            case '1D1D'
                y0 = [p20; dp0; c0; dpp0; cp0]; % pack initial conditions
        end
        
        tBeg = tEnd; % start time of simulation [s]
        tEnd = tBeg + 2*t0; % end time of simulation [s]

        % Integrate in time -------------------------------------------------------
        [tt_t,yy] = ode15s(@(t,y) rhsCV(t,y,par,parIn), tBeg:tStep:tEnd, y0, options);
        
        phi2_t = yy(:,1:2*N+Ns); % potential in electrolyte (phi2, algebraic variable)
        dphi_t = yy(:,2*N+Ns+1:4*N+Ns); % potential difference (phi1-phi2)
        c_t = yy(:,4*N+Ns+1:6*N+2*Ns); % concentration
        
        % extract phi1 (potential in the solid phase)
        phi1_t = NaN*phi2_t;
        phi1_t(:,1:N) = dphi_t(:,1:N) + phi2_t(:,1:N);
        phi1_t(:,N+Ns+1:2*N+Ns) = dphi_t(:,N+1:2*N) + phi2_t(:,N+Ns+1:2*N+Ns);
        
        U_t = NaN*tt_t; % initialize
        for i=1:length(tt_t)
            U_t(i,1) = UMax - nuVec(k)*abs(mod(tt_t(i),2*t0)-t0); % compute voltage
        end
        
        j_t = 2.0*sigmaS*phi1_t(:,1)/h; % current density at the left current collector
        
        qs = trapz(U_t,j_t*S/2/nuVec(k)); % compute charge
        Cs_n(n) = qs/(UMax - UMin); % compute capacitance

        phiLo = min(dphi_t(:,1:2*N),[],1); % minimum potential during simulation
        phiHi = max(dphi_t(:,1:2*N),[],1); % maximum potential during simulation
        
        sumL = 1/(N-1)*trapz(phiHi(1:N)-phiLo(1:N)); % summation
        sumR = 1/(N-1)*trapz(phiHi(N+1:2*N)-phiLo(N+1:2*N)); % summation

        lambda(n) = (sumL+sumR)/(UMax-UMin); % penetration depth

        disp(['C = ' num2str(Cs_n(n)) ' F'])
        disp(['lambda = ' num2str(lambda(n))])
        
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

                phi2_p((n-1)*n_t+1 : n*n_t,:) = phi2_p_t(:,:);
                cp((n-1)*n_t+1 : n*n_t,:) = cp_t(:,:);
        end
        
        % store this cycle in the overall vectors
        tt((n-1)*n_t+1 : n*n_t) = tt_t(:); 
        phi2((n-1)*n_t+1 : n*n_t,:) = phi2_t(:,:);
        phi1((n-1)*n_t+1 : n*n_t,:) = phi1_t(:,:);
        c((n-1)*n_t+1 : n*n_t,:) = c_t(:,:);
        jj((n-1)*n_t+1 : n*n_t) = j_t(:);
        UU((n-1)*n_t+1 : n*n_t) = U_t(:);
        cMin((n-1)*n_t+1 : n*n_t) = cMin_t(:);
        cMax((n-1)*n_t+1 : n*n_t) = cMax_t(:);
        
        % conditional break if steady-state capacitance is reached
        if n~=1
            if abs(Cs_n(n) - Cs_n(n-1)) < thresh*Cmax;
                % shorten vectors
                tt = tt(1 : n*n_t);
                phi2 = phi2(1 : n*n_t,:);
                phi1 = phi1(1 : n*n_t,:);
                c = c(1 : n*n_t,:);
                jj = jj(1 : n*n_t);
                UU = UU(1 : n*n_t);
                cMin = cMin(1 : n*n_t);
                cMax = cMax(1 : n*n_t);
                Cs_n = Cs_n(1:n);
                lambda = lambda(1:n);
                break
            end
        end

    end
    
    Cs(k) = Cs_n(end);
    
    % save fields to exported variables (structure res)
    CV.nu{k}.tt = tt;
    CV.nu{k}.phi2 = phi2;
    CV.nu{k}.phi1 = phi1;
    CV.nu{k}.c = c;
    CV.nu{k}.jj = jj;
    CV.nu{k}.UU = UU;
    CV.nu{k}.Cs = Cs(k);
    CV.nu{k}.lambda = lambda(end);
    CV.nu{k}.cMin = cMin;
    CV.nu{k}.cMax = cMax;
    
    switch dim
        case '1D1D'
            CV.nu{k}.phi2_p = phi2_p;
            CV.nu{k}.cp = cp;
    end
    
end

% compute time constant
td = NaN; % reset
Cs_max = max(Cs);
for k=1:N_nu
    if Cs(k) <= Cs_max/2
        % interpolate knee scan rate
        nuK = nuVec(k-1) + (nuVec(k) - nuVec(k-1))* ...
            (Cs_max/2 - Cs(k-1))/(Cs(k) - Cs(k-1)); % lever rule
        td = 2*(UMax - UMin)/nuK; % time constant
        break
    end
end

CV.td = td;






















