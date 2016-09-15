% This file is part of the program SupPer for the modelling of
% dynamic behavior and performance of Supercapacitors
% (Electrochemical Double Layer Capacitors)
% 
% Author: Martin Kroupa
% Oct 2015 - Mar 2016
% 
% Developed jointly at Imperial College London and University of Chemistry
% and Technology Prague

%% function for plotting capacity as a function of frequency
function figureCapFreq(res, par)

%% Parameters
v2struct(par.fig);
v2struct(par.gen);

%% Plotted variables

if isfield(res,'CV')
    v2struct(par.CV);
    fVec_CV = nuVec/(2*(UMax-UMin));
    Cs_CV = NaN*zeros(N_nu,1);
    lambda_CV = NaN*zeros(N_nu,1);
    for k=1:N_nu
        Cs_CV(k) = res.CV.nu{k}.Cs;
        lambda_CV(k) = res.CV.nu{k}.lambda;
    end
else
    fVec_CV = NaN;
    Cs_CV = NaN;
end

if isfield(res,'EIS')
    v2struct(par.EIS);
    fVec_EIS = fVec;
    Cs_EIS = NaN*zeros(N_f,1);
    lambda_EIS = NaN*zeros(N_f,1);
    for k=1:N_f
        Cs_EIS(k) = res.EIS.f{k}.Cs;
        lambda_EIS(k) = res.EIS.f{k}.lambda;
    end
else
    fVec_EIS = NaN;
    Cs_EIS = NaN;
end

if isfield(res,'GC')
    v2struct(par.GC);
    fVec_GC = NaN*zeros(N_I,1);
    Cs_GC = NaN*zeros(N_I,1);
    lambda_GC = NaN*zeros(N_I,1);
    for k=1:N_I
        fVec_GC(k) = res.GC.I{k}.f;
        Cs_GC(k) = res.GC.I{k}.Cs;
        lambda_GC(k) = res.GC.I{k}.lambda;
    end
else
    fVec_GC = NaN;
    Cs_GC = NaN;
end

%% plot capacity versus frequency

xPos = xSize;
yPos = 1;

figure1 = figure;
set(figure1, 'Units','centimeters')
set(figure1, 'Position', [xPos yPos xSize ySize])
pos = get(figure1,'Position');
set(figure1, 'PaperPositionMode' , 'Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])

plot(fVec_CV,Cs_CV,'bo-','MarkerSize', 10,'LineWidth', line2);
hold on
plot(fVec_EIS,Cs_EIS,'rx--','MarkerSize', 10,'LineWidth', line2);
plot(fVec_GC,Cs_GC,'kd-.','MarkerSize', 10,'LineWidth', line2);

grid on
% set(gca,'YScale','log')
set(gca,'XScale','log')
% 
% set(gca,'ylim',[0 8e-5])
set(gca,'xlim',[1e-4 1e2])
set(gca,'xTick',[1e-4 1e-3 1e-2 1e-1 1 1e1 1e2])

xlabel('Frequency ({\itf} / s^{-1})','FontSize', fontSize) 
ylabel('Real capacitance ({\itC}^{Re} / F)','FontSize', fontSize) 
set(gca,'FontSize', fontSize);
