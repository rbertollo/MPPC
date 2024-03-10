%% Observability and Identifiability Tests
clc
clear all
close all

addpath('/Users/marcodelloro/Downloads/casadi-3.6.3-osx64-matlab2018b')
import casadi.*;
opti = casadi.Opti();

load('/Users/marcodelloro/Desktop/Thesis/MSc-Thesis-TUe/Data_Collection/Australian Dataset/SIDTTHE_AUS.mat');
load('/Users/marcodelloro/Desktop/Thesis/MSc-Thesis-TUe/Data_Collection/Australian Dataset/DataStructAUS.mat');
set(0,'DefaultFigureWindowStyle','docked');

%%  Data Loading and Initialization

Npop = 25766605; % Total Population of Australia as of 31 December 2021
N = 399;       % MPC horizon length
N_mhe = 21;    % Estimation horizon (3 weeks)
date = DataStructAUS.Raw.Date;

I_data = SIDTTHE_AUS{1,1} ./ Npop;     
D_data = SIDTTHE_AUS{2,1} ./ Npop;     
T_data = (SIDTTHE_AUS{3,1} + SIDTTHE_AUS{4,1}) ./ Npop;
H_data = SIDTTHE_AUS{7,1} ./ Npop;
H_dataAug = SIDTTHE_AUS{6,1}' ./ Npop;
E_data = SIDTTHE_AUS{5,1}  ./ Npop;
S_data = ones(length(I_data),1) - (I_data + D_data + T_data + H_dataAug + E_data );

%% Define the set of differential equations and relative syms variables

syms S I D T H E % Variables
syms alpha gamma delta lambda lambda2 sigma tau % Coefs

eqns = [    -S * (alpha * I);...
             S * (alpha * I) - (gamma + lambda) * I;...
             I * gamma - D * (lambda + delta);...
             delta * D - (tau + sigma)*T;...
             lambda * D + T * sigma + lambda * I    ];


vars = [S I D T H alpha gamma delta sigma tau];
coefs = [lambda];

% Passing the vector of equations 

[f1,f2,f_ODE] = getdynamics(eqns,vars,coefs);

% Write "C" measurements matrix
C = [ diag([1 1 1 1 1]), zeros(5,5)];

J = jacobian(f1,vars);

Ob = obsv_sym(J,C);
if rank(Ob) == size(J,1)
    disp("OBSERVABLE, Full Rank")
else
    fprintf('Non observable, rank: %d', rank(Ob));
end

%% Fitting in order to find the value of the unknown parameter

N = 69; % timeframe of the fit

DataArray = { S_data(1:N)' ; I_data(1:N)' ; D_data(1:N)' ; T_data(1:N)' ; H_data(1:N)' };

data = cell2mat(DataArray);
date = DataStructAUS.Raw.Date;

S = opti.variable(1,N); 
I = opti.variable(1,N);
D = opti.variable(1,N);
T = opti.variable(1,N);
H = opti.variable(1,N);
alpha = opti.variable(1,N);
gamma = opti.variable(1,N);
delta = opti.variable(1,N);
sigma = opti.variable(1,N);
tau = opti.variable(1,N);
lambda = opti.variable(1,N);

x = [S; I; D; T; H; alpha; gamma; delta; sigma; tau];
par = [lambda];

% Constraints related to the coefficients

opti.subject_to(lambda(1:end) >= 0);   % bound on lambda value
opti.subject_to(lambda(1:end) <= 0.3);

for jj = 1:N-1
        opti.subject_to( lambda(1,jj + 1) == lambda(1,jj) )
end

% Simulation with RK45

dt = 1;

for k = 1:N-1 
   
   k1 = f2(x(:,k),         par(:,k));
   k2 = f2(x(:,k)+dt/2*k1, par(:,k));
   k3 = f2(x(:,k)+dt/2*k2, par(:,k));
   k4 = f2(x(:,k)+dt*k3,   par(:,k));
   x_next = x(:,k) + dt/6*(k1+2*k2+2*k3+k4);
   opti.subject_to(x(:,k+1) == x_next); % close the gaps
end

% Optimization 

data_obj = horzcat(data(:)); 
x_obj = [x(1, :), x(2, :), x(3, :), x(4, :), x(5, :)];

maxData =( ones(N,1)*max(data,[],2)' )';
maxData = horzcat(maxData(:));

coefs_matrix_obj = [ sum( (diff(lambda)./0.1).^2 )];

obj = sum(((data_obj - x_obj')./maxData).^2)*0.01 + 1e-3*coefs_matrix_obj;
opti.minimize(obj);
p_opts = struct('expand', false);
s_opts = struct('max_iter',1e4, 'tol',1e-4,'constr_viol_tol',1e-4,'compl_inf_tol',1e-4,'linear_solver','MA97'); %'MA27' 'MA57' 'MA77' 'MA86' 'MA97'
opti.solver('ipopt',p_opts, s_opts); % set numerical backend
sol = opti.solve();   % actual solver

column_names = {'lambda'};
opti_lambda = array2table(opti.debug.value(lambda)', 'VariableNames', column_names);

% Plot of the fitted coefficient
figure(1)
plot(date(1:N),opti_lambda.lambda, LineWidth=1.5)
ylabel('Coefficient Value','Interpreter','latex')
title('\textbf{$\lambda$ coefficient}','Interpreter','latex')
grid on
xlim([date(1), date(N)])
ylim([0, max(opti_lambda.lambda)*1.05])
set(gca, 'TickLabelInterpreter', 'Latex')