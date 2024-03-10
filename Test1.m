%% Observability and Identifiability Tests
clear all
clc

% Define the set of differential equations and relative syms variables

syms S I D T1 T2 H E % Variables
syms alpha gamma delta1 delta2 lambda sigma1 sigma2 tau1 tau2 % Coefs


eqns = [    -S * (alpha * I);...
             S * (alpha * I) - (gamma + lambda) * I;...
             I * gamma - D * (lambda + delta1 + delta2);...
             delta1 * D - (sigma1 + tau1) * T1;...
             delta2 * D - (sigma2 + tau2) * T2;...
             (I + D) * lambda + T1 * sigma1 + T2 * sigma2 ];

% eqns = [    -S * (alpha * I) + lambda * I;...
%              S * (alpha * I) - (gamma + lambda) * I;...
%              I * gamma - D * (lambda + delta1 + delta2);...
%              delta1 * D - (sigma1 + tau1) * T1;...
%              delta2 * D - (sigma2 + tau2) * T2;...
%              D * lambda + T1 * sigma1 + T2 * sigma2 ];

vars = [S I D T1 T2 H alpha gamma delta1 delta2 sigma1 tau1];
coefs = [sigma2 tau2 lambda];

% Passing the vector of equations 

[f1,f2] = getdynamics(eqns,vars,coefs);

% Write "C" measurements matrix
C = [eye(6), zeros(6,6)];

J = jacobian(f1,vars);

Ob = obsv_sym(J,C);
if rank(Ob) == size(J,1)
    disp("OBSERVABLE, Full Rank")
else
    fprintf('Non observable, rank: %d', rank(Ob));
end

%% Identifiability test
% Load of the data 

load('/Users/marcodelloro/Desktop/Thesis/MSc-Thesis-TUe/Data_Collection/Italian Dataset/SIDTTHE_data_DEF.mat');
date = SIDTTHE_data{1,1}.date;
Npop = 59240329; % Total Population of Italy

I_data = SIDTTHE_data{1,1}.data / Npop;     
D_data = SIDTTHE_data{2,1}.data / Npop;     
T1_data = SIDTTHE_data{3,1}.data / Npop;
T2_data = SIDTTHE_data{4,1}.data / Npop;
H_dataAug = SIDTTHE_data{6,1}.data / Npop;
H_data = SIDTTHE_data{7,1}.data / Npop;
E_data = SIDTTHE_data{5,1}.data  / Npop;
S_data = ones(length(I_data),1)' - (I_data + D_data + T1_data + T2_data + H_data + E_data );

% Initial values for parameters both in dynamics and estimation
alpha_0 = 0.513;
gamma_0 = 0.24;
delta1_0 = 0.01;
delta2_0 = 0.005;
sigma1_0 = 0.03;
sigma2_0 = 0.065;
tau1_0 = 0.015;
tau2_0 = 0.005;
lambda_0 = 0.0596;

x0 = [S_data(1); I_data(1); D_data(1); T1_data(1); T2_data(1); H_data(1); alpha_0; gamma_0; delta1_0; delta2_0; sigma1_0; tau1_0;];
par0 = [lambda_0; sigma2_0; tau2_0];
t = 1:399;
data = [S_data', I_data', D_data', T1_data', T2_data', H_data'];

[opt, DeltaV] = identif_test(f2, C, t', data, x0, par0);
DeltaV