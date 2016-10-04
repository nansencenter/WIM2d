function [E,g,rho_wtr,rho_ice,nu]=NDphyspram_v2();
%% THIS PROGRAM STORES ALL THE PHYSICAL PARAMETERS NEEDED FOR THE PROBLEM.
%% CALL: [E,g,rho_wtr,rho_ice,nu]=NDphyspram;

E        =5e9;%% Pa
g        =9.81;%% m/s^2
rho_wtr  =1025;%% kg/m^3
rho_ice  =922.5;%% kg/m^3
nu       =.3;
