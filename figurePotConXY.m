% This file is part of the program SupPer for the modelling of
% dynamic behavior and performance of Supercapacitors
% (Electrochemical Double Layer Capacitors)
% 
% Author: Martin Kroupa
% Oct 2015 - Mar 2016
% 
% Developed jointly at Imperial College London and University of Chemistry
% and Technology Prague

%% function for plotting potential as a function of time
function figurePotConXY(res, par, k, tEnd)

%% Parameters
v2struct(par.fig);
v2struct(par.num);

method = par.spec.method;

%% Plotted variables

rx = (1:N)*h - 0.5*h;
rx(N+1) = rx(N) + 0.5*h + 0.5*hs;
for i=2:Ns
rx(N+i) = rx(N+i-1) + hs;
end
rx(N+Ns+1) = rx(N+Ns) + 0.5*hs + 0.5*h;
for i=N+Ns+2:2*N+Ns
rx(i) = rx(i-1) + h;
end
rx = rx*1e6; % convert to micrometers

ry = (1:Np)*hp - 0.5*hp;
ry = ry*1e6; % convert to micrometers

switch method
    case 'CV'
        tt = res.nu{k}.tt;
        phi2_p = res.nu{k}.phi2_p;
        cp = res.nu{k}.cp;
    case 'EIS'
        tt = res.f{k}.tt;
        phi2_p = real(res.f{k}.phi2_p);
        cp = real(res.f{k}.cp);
    case 'GC'
        tt = res.I{k}.tt;
        phi2_p = res.I{k}.phi2_p;
        cp = res.I{k}.cp;
    case 'CCC'
        tt = res.tt;
        phi2_p = res.phi2_p;
        cp = res.cp;
end

tEndIdx = find(tt<=tEnd);
tEndIdx = tEndIdx(end);

%% plot concentration profile evolution

xPos = 0;
yPos = 17;

figure1 = figure;
set(figure1, 'Units','centimeters')
set(figure1, 'Position', [xPos yPos xSize ySize])
pos = get(figure1,'Position');
set(figure1, 'PaperPositionMode' , 'Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])

cp_plot = NaN*zeros(Np,2*N+Ns);

for i=1:N
    for j=1:Np
       cp_plot(j,i) = cp(tEndIdx,(i-1)*Np+j); 
    end
end

for i=N+1:2*N
    for j=1:Np
       cp_plot(j,i+Ns) = cp(tEndIdx,(i-1)*Np+j); 
    end
end

surf(rx,ry,cp_plot);

% set(gca,'YScale','log')
% set(gca,'XScale','log')
% 
% set(gca,'ylim',[0 0.8])
% set(gca,'xlim',[0 125])

xlabel('x-coordinate ({\itx} / \mum)','FontSize', fontSize)
ylabel('y-coordinate ({\ity} / \mum)','FontSize', fontSize) 
zlabel('concentration in meso-pores ({\it\phi}_l / V)','FontSize', fontSize) 
set(gca,'FontSize', fontSize);

%% plot concentration profile evolution

xPos = xSize;
yPos = 17;

figure1 = figure;
set(figure1, 'Units','centimeters')
set(figure1, 'Position', [xPos yPos xSize ySize])
pos = get(figure1,'Position');
set(figure1, 'PaperPositionMode' , 'Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])

phi2_p_plot = NaN*zeros(Np,2*N+Ns);

for i=1:N
    for j=1:Np
       phi2_p_plot(j,i) = phi2_p(tEndIdx,(i-1)*Np+j); 
    end
end

for i=N+1:2*N
    for j=1:Np
       phi2_p_plot(j,i+Ns) = phi2_p(tEndIdx,(i-1)*Np+j); 
    end
end

surf(rx,ry,phi2_p_plot);

% set(gca,'YScale','log')
% set(gca,'XScale','log')
% 
% set(gca,'ylim',[0 0.8])
% set(gca,'xlim',[0 125])

xlabel('x-coordinate ({\itx} / \mum)','FontSize', fontSize)
ylabel('y-coordinate ({\ity} / \mum)','FontSize', fontSize) 
zlabel('potential in meso-pores ({\it\phi}_l / V)','FontSize', fontSize) 
set(gca,'FontSize', fontSize);

end











