% This file is part of the program SupPer for the modelling of
% dynamic behavior and performance of Supercapacitors
% (Electrochemical Double Layer Capacitors)
% 
% Author: Martin Kroupa
% Oct 2015 - Mar 2016
% 
% Developed jointly at Imperial College London and University of Chemistry
% and Technology Prague

%% function for plotting voltage and Current Characteristics
function figureCurVolTime(res, par, k)

%% Parameters
v2struct(par.fig);
v2struct(par.gen);
v2struct(par.spec); % unpack specification parameters

%% Plotted variables
switch method
    case 'CV'
        tt = res.nu{k}.tt;
        UU = res.nu{k}.UU;
        jj = res.nu{k}.jj;
        cMin = res.nu{k}.cMin;
        cMax = res.nu{k}.cMax;
    case 'EIS'
        tt = res.f{k}.tt;
        UU = real(res.f{k}.UU);
        jj = real(res.f{k}.jj);
        cMin = real(res.f{k}.cMin);
        cMax = real(res.f{k}.cMax);
    case 'GC'
        tt = res.I{k}.tt;
        UU = res.I{k}.UU;
        jj = res.I{k}.jj;
        cMin = res.I{k}.cMin;
        cMax = res.I{k}.cMax;
    case 'CCC'
        tt = res.tt;
        UU = res.UU;
        jj = res.jj;
        cMin = res.cMin;
        cMax = res.cMax;
end

%% plot voltage versus time

xPos = 0;
yPos = 16;

figure1 = figure;
set(figure1, 'Units','centimeters')
set(figure1, 'Position', [xPos yPos xSize ySize])
pos = get(figure1,'Position');
set(figure1, 'PaperPositionMode' , 'Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])

plot(tt,UU,'k-','MarkerSize', 10,'LineWidth', line2);

grid on
% set(gca,'YScale','log')
% set(gca,'XScale','log')
% 
% set(gca,'ylim',[1.4 2.8])
% set(gca,'yTick',[1.4 1.6 1.8 2.0 2.2 2.4 2.6 2.8])
% set(gca,'xlim',[0 800])

xlabel('Time ({\itt} / s)','FontSize', fontSize)
ylabel('Voltage ({\itU} / V)','FontSize', fontSize) 
set(gca,'FontSize', fontSize);


%% plot current versus time

xPos = xSize;
yPos = 16;

figure1 = figure;
set(figure1, 'Units','centimeters')
set(figure1, 'Position', [xPos yPos xSize ySize])
pos = get(figure1,'Position');
set(figure1, 'PaperPositionMode' , 'Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])

plot(tt,jj*S,'bo','MarkerSize', 10,'LineWidth', line2);

grid on
% set(gca,'YScale','log')
% set(gca,'XScale','log')
% 
% set(gca,'ylim',[0 0.8])
% set(gca,'xlim',[0 125])

xlabel('Time ({\itt} / s)','FontSize', fontSize)
ylabel('Current ({\itI} / A)','FontSize', fontSize) 
set(gca,'FontSize', fontSize);


%% plot current versus potential

xPos = 2*xSize;
yPos = 16;

figure1 = figure;
set(figure1, 'Units','centimeters')
set(figure1, 'Position', [xPos yPos xSize ySize])
pos = get(figure1,'Position');
set(figure1, 'PaperPositionMode' , 'Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])

line1 = 1;
plot(UU,jj*S,'rx','MarkerSize', 5,'LineWidth', line1);

grid on
% set(gca,'YScale','log')
% set(gca,'XScale','log')
% 
% set(gca,'ylim',[0 0.8])
% set(gca,'xlim',[0 125])

xlabel('Voltage ({\itU} / V)','FontSize', fontSize)
ylabel('Current ({\itI} / A)','FontSize', fontSize) 
set(gca,'FontSize', fontSize);

%% plot concentration extremes versus time

xPos = 0;
yPos = 1;

figure1 = figure;
set(figure1, 'Units','centimeters')
set(figure1, 'Position', [xPos yPos xSize ySize])
pos = get(figure1,'Position');
set(figure1, 'PaperPositionMode' , 'Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])

plot(tt,cMin,'k-','MarkerSize', 10,'LineWidth', line2);
hold on
plot(tt,cMax,'r--','MarkerSize', 10,'LineWidth', line2);

grid on
% set(gca,'YScale','log')
% set(gca,'XScale','log')
% 
% set(gca,'ylim',[9.99 10.01])
% set(gca,'xlim',[0 125])

xlabel('Time ({\itt} / s)','FontSize', fontSize)
ylabel('Concentration ({\itc} / mol\cdotm^{-3})','FontSize', fontSize) 
set(gca,'FontSize', fontSize);

legend('minimum'...
    ,'maximum'...
    ,'Location','SouthEast');


