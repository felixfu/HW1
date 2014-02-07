clear all
clc

%============Inputs============

%**System Paremeters**
TDP = 125.0;           % Thermal Design Power, W
NumFins = 5;           % Number of heat sink fins
L = 8/100;             % Total length of domain (CPU + heat sink), m
L_CPU = 3/100;         % Length of CPU region, m
thick = 0.3/100;       % Thickness of heat sink fin, m

%**Material Parameters**
k_Metal = 160.;               % Thermal conductivity of Metal, W/(m*K)
density_Metal = 2700.;        % Density of metal, kg/m^3
specific_heat_Metal = 895.0;  % Specific heat of metal, J/(kg*K)
alpha_Metal = k_Metal/(density_Metal*specific_heat_Metal);   % Thermal diffusivity of Metal, m^2/s
%====NOTE: k_CPU is the uncertain input variable
%===========================================================
k_CPU =0.40;                  % Thermal conductivity of CPU, W/(m*K)
%===========================================================
density_CPU = 2000.;          % Density of CPU, kg/m^3
specific_heat_CPU = 750.;     % Specific heat of CPU, J/(kg*K)
alpha_CPU = k_CPU/(density_CPU*specific_heat_CPU); % Thermal diffusivity of CPU, m^2/s

%**Model Parameters (averaged in flow direction)**
hbar = 15.0;           % Averaged convective cooling coeff., W/(m^2*K)
Tbar = 330.0;          % Averaged air temperature, K

%**Discretization Parameters**
imax = 33;             % Number of nodes in x (normal to flow direction) 
itermax = 1e7;         % Maximum allowable number of iterations
convtol = 1.e-8;       % Iterative convergence tolerance (relative to fifth iteration)

%**Solution Variable**
T(1:imax) = 350.0;     % Temperature of heat sink (initialized to 350 K), K

[Tbase, T, x, L2conv, history] = heatcondsolve(T,TDP,NumFins,L,thick,L_CPU,k_Metal,alpha_Metal,k_CPU,alpha_CPU,hbar,Tbar,imax,itermax,convtol);

Tbase-273.15

% NOTE: you may want to suppress figure output for nondeterministic simulations
figure(1)
subplot(2,1,1)
semilogy(history(:,1),history(:,2))
title('Iterative convergence')

subplot(2,1,2)
plot(x,T-273.15,'b-o')
title('Temperature')





