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

vars = [S I D T1 T2 H alpha gamma delta1 delta2 sigma1 sigma2 tau1 tau2];
coefs = [lambda];

% Passing the vector of equations 

[f1,f2] = getdynamics(eqns,vars,coefs);

% Write "C" measurements matrix
C = [eye(6), zeros(6,8)];

J = jacobian(f1,vars);

Ob = obsv_sym(J,C);
if rank(Ob) == size(J,1)
    disp("OBSERVABLE, Full Rank")
else
    fprintf('Non observable, rank: %d', rank(Ob));
end
