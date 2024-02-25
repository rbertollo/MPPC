%% Observability and Identifiability Tests of SIDARTHE-V model
clear all
close all
clc

% Define the set of differential equations and relative syms variables

sym S I D A R T H E 
sym alpha beta gamma delta epsilon zeta eta theta mu nu kappa lambda sigma tau1 tau2 phi psi rho

% Defining the equations
eqns = [   -S*(alpha*I + beta*D + gamma*A + delta*R);...
            S*(alpha*I + beta*D + gamma*A + delta*R) - (epsilon + zeta + lambda)*I;...
            epsilon*I - (eta + rho)*D;...
            eta*I - (theta + mu + kappa)*A;...
            eta*D + theta*A - (nu + psi + tau1)*R;...
            mu*A + nu*R - (sigma + tau2)*T;...
            lambda*I + rho*D + kappa*A + psi*R + sigma*T;...
            tau1*R + tau2*T                         ];

vars = [S I D A R T H E alpha beta gamma delta epsilon sigma tau1 tau2 phi psi rho];
coefs = [ zeta eta theta mu nu kappa lambda];

% Passing the vector of equations 

[f1,~,~] = getdynamics(eqns,vars,coefs);

% Write "C" measurements matrix
C = [eye(8), zeros(8,11)];

J = jacobian(f1,vars);

Ob = obsv_sym(J,C);
if rank(Ob) == size(J,1)
    disp("OBSERVABLE, Full Rank")
else
    fprintf('Non observable, rank: %d', rank(Ob));
end
