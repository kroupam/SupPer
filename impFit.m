function F = impFit(cc,tt,j,f)
% This file is part of the program SupPer for the modelling of
% dynamic behavior and performance of Supercapacitors
% (Electrochemical Double Layer Capacitors)
% 
% Author: Martin Kroupa
% Oct 2015 - Mar 2016
% 
% Developed jointly at Imperial College London and University of Chemistry
% and Technology Prague

% Function for the fitting of impedance from U vs. t and I vs. t data

% cc is a two-element vector whose components are the real
% and imaginary parts of the parameter respectively

j0 = cc(1);
phi_s = cc(2);

errors = j0*exp(1i*(2*pi*f*tt+phi_s)) - j;

F = [real(errors(:)); imag(errors(:))];