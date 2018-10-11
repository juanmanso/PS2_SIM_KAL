%
% Archivo de configuraci√≥n del main para el graph.m
%

clear all, close all;

%EsMatlab = 1;
EsMatlab = sum(sum(ismember(struct2cell(ver),'MATLAB'))); % Chequeo autom·tico de MATLAB
if(EsMatlab == 0)
    graphics_toolkit('gnuplot');
end

format long;

% Agrego las funcinoes definidas en la carpeta "funciones"
addpath('./Funciones')
addpath('./Material')

myGreen=[0 0.5 0];

