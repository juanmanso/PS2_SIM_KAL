%
% Archivo de configuración del main para el graph.m
%

clear all, close all;

EsMatlab = 1;
if(EsMatlab == 0)
    graphics_toolkit("gnuplot");
end

format long;

% Agrego las funcinoes definidas en la carpeta "funciones"
addpath('./Funciones')
addpath('./Material')

myGreen=[0 0.5 0];

