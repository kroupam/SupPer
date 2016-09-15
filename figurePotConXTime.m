% This file is part of the program SupPer for the modelling of
% dynamic behavior and performance of Supercapacitors
% (Electrochemical Double Layer Capacitors)
% 
% Author: Martin Kroupa
% Oct 2015 - Mar 2016
% 
% Developed jointly at Imperial College London and University of Chemistry
% and Technology Prague

%% function for plotting potential as a function of space and time
function figurePotConXTime(res, par, k, tBeg, tEnd)

%% Parameters
v2struct(par.fig);
v2struct(par.num);

method = par.spec.method;

plotStep = 15; % plot with coarser resolution

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

switch method
    case 'CV'
        tt = res.nu{k}.tt;
        phi2 = res.nu{k}.phi2;
        phi1 = res.nu{k}.phi1;
        c = res.nu{k}.c;
    case 'EIS'
        tt = res.f{k}.tt;
        phi2 = real(res.f{k}.phi2);
        phi1 = real(res.f{k}.phi1);
        c = real(res.f{k}.c);
    case 'GC'
        tt = res.I{k}.tt;
        phi2 = res.I{k}.phi2;
        phi1 = res.I{k}.phi1;
        c = res.I{k}.c;
    case 'CPC'
        tt = res.P{k}.tt;
        phi2 = res.P{k}.phi2;
        phi1 = res.P{k}.phi1;
        c = res.P{k}.c;
    case 'CCC'
        tt = res.tt;
        phi2 = res.phi2;
        phi1 = res.phi1;
        c = res.c;
end

tBegIdx = find(tt<=tBeg);
tBegIdx = tBegIdx(end);
tEndIdx = find(tt<=tEnd);
tEndIdx = tEndIdx(end);

%% plot potential in electrolyte profile evolution

xPos = 0;
yPos = 2;

figure1 = figure;
set(figure1, 'Units','centimeters')
set(figure1, 'Position', [xPos yPos xSize ySize])
pos = get(figure1,'Position');
set(figure1, 'PaperPositionMode' , 'Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])

p1 = plot(rx,phi2(tBegIdx,:),'k-','MarkerSize', 10,'LineWidth', line2);
hold on
for i=tBegIdx+plotStep:plotStep:tEndIdx-1
    plot(rx,phi2(i,:),'r-')
end
p_end = plot(rx,phi2(tEndIdx,:),'bo-','MarkerSize', 10,'LineWidth', line2);
hold off

grid on
% set(gca,'YScale','log')
% set(gca,'XScale','log')
% 
set(gca,'ylim',[0.3 1.2])
set(gca,'xlim',[0 125])
set(gca,'xTick',[0 25 50 75 100 125])

xlabel('{\itx}-coordinate ({\itx} / \mum)','FontSize', fontSize)
ylabel('Potential in electrolyte (\Phi_2 / V)','FontSize', fontSize) 
set(gca,'FontSize', fontSize);

legend([ ...
    p1 ...
    p_end ...
    ],{ ...
    ['t = ' num2str(tBeg) ' s'] ...
    ,['t = ' num2str(tEnd) ' s']} ...
    ,'Location','SouthEast');

%% plot solid-phase potential profile evolution

xPos = xSize;
yPos = 2;

figure1 = figure;
set(figure1, 'Units','centimeters')
set(figure1, 'Position', [xPos yPos xSize ySize])
pos = get(figure1,'Position');
set(figure1, 'PaperPositionMode' , 'Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])

p1 = plot(rx,phi1(tBegIdx,:),'k-','MarkerSize', 10,'LineWidth', line2);
hold on
for i=tBegIdx+plotStep:plotStep:tEndIdx-1
    plot(rx,phi1(i,:),'r-')
end
p_end = plot(rx,phi1(tEndIdx,:),'bo-','MarkerSize', 10,'LineWidth', line2);
hold off

grid on
% set(gca,'YScale','log')
% set(gca,'XScale','log')
% 
set(gca,'ylim',[-0.5 2.5])
set(gca,'xlim',[0 125])
set(gca,'xTick',[0 25 50 75 100 125])

xlabel('{\itx}-coordinate ({\itx} / \mum)','FontSize', fontSize)
ylabel('Solid-phase potential (\Phi_1 / V)','FontSize', fontSize) 
set(gca,'FontSize', fontSize);

legend([ ...
    p1 ...
    p_end ...
    ],{ ...
    ['t = ' num2str(tBeg) ' s'] ...
    ,['t = ' num2str(tEnd) ' s']} ...
    ,'Location','SouthEast');

%% plot concentration profile evolution

xPos = 2*xSize;
yPos = 2;

figure1 = figure;
set(figure1, 'Units','centimeters')
set(figure1, 'Position', [xPos yPos xSize ySize])
pos = get(figure1,'Position');
set(figure1, 'PaperPositionMode' , 'Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])

c1 = plot(rx,c(tBegIdx,:),'k-','MarkerSize', 10,'LineWidth', line2);
hold on
for i=tBegIdx+plotStep:plotStep:tEndIdx-1
    plot(rx,c(i,:),'r-')
end
c_end = plot(rx,c(tEndIdx,:),'bo-','MarkerSize', 10,'LineWidth', line2);
hold off

grid on
% set(gca,'YScale','log')
% set(gca,'XScale','log')
% 
% set(gca,'ylim',[800 1000])
set(gca,'xlim',[0 125])
set(gca,'xTick',[0 25 50 75 100 125])

xlabel('{\itx}-coordinate ({\itx} / \mum)','FontSize', fontSize)
ylabel('Concentration of salt ({\itc} / mol\cdotm^{-3})','FontSize', fontSize) 
set(gca,'FontSize', fontSize);

legend([ ...
    c1 ...
    c_end ...
    ],{ ...
    ['t = ' num2str(tBeg) ' s'] ...
    ,['t = ' num2str(tEnd) ' s']} ...
    ,'Location','SouthEast');


end











