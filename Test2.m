%% Observability and Identifiability Tests
clear all
clc


addpath('/Users/marcodelloro/Downloads/casadi-3.6.3-osx64-matlab2018b')
import casadi.*;
opti = casadi.Opti();

% Define the set of differential equations and relative syms variables

syms S I D T H E % Variables
syms alpha gamma delta lambda sigma tau % Coefs

eqns = [    -S * (alpha * I);...
             S * (alpha * I) - (gamma+lambda) * I;...
             I * gamma - D * (lambda + delta);...
             delta * D - (tau + sigma)*T
             lambda * (D+I) + T * sigma ];

vars = [S I D T H alpha gamma delta sigma tau];
coefs = [lambda];

% Passing the vector of equations 

[f1,f2,f_ODE] = getdynamics(eqns,vars,coefs);

% Write "C" measurements matrix
C = [eye(5), zeros(5,5)];

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

T_data = T1_data + T2_data;

% Initial values for parameters both in dynamics and estimation
alpha_0 = 0.26;
gamma_0 = 0.24;
delta_0 = 0.01;
sigma_0 = 0.03;
tau_0 = 0.01;
lambda_0 = 0.0596;

x0 = [S_data(1); I_data(1); D_data(1); T_data(1); H_data(1); alpha_0; gamma_0; delta_0; sigma_0; tau_0;];
par0 = [lambda_0];
t = 1:399;
data = [S_data', I_data', D_data', T_data', H_data'];

[opt, DeltaV] = identif_test(f2, C, t', data, x0, par0);

%% Fitting in order to find the value of the unknown parameter

N = 69; % timeframe of the fit

DataArray = {S_data(1,1:N) ; I_data(1,1:N) ; D_data(1,1:N) ; T_data(1,1:N) ; H_data(1,1:N)};

data = cell2mat(DataArray);
date = SIDTTHE_data{1,1}.date;

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
opti.subject_to(lambda(1:end) <= 0.2);

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

obj = sum(((data_obj - x_obj')./maxData).^2)*0.01;
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

%% EKF computation 

% Raw data loading for std
dataset = readtable('Italy_complete_dataset.xlsx');
dataset = dataset(190:588,:);
dataset.data = datetime(dataset.data, 'InputFormat', 'dd-MMM-yyyy HH:mm:ss');

% Normalisation of data
I_data = I_data./max(I_data);
D_data = D_data./max(D_data);
T_data = T_data./max(T_data);
H_data = H_data./max(H_data);
E_data = E_data./max(E_data);
S_data = S_data./max(S_data);

rawT = dataset.ricoverati + dataset.terapia_intensiva/Npop;
rawD = dataset.isolamento_domiciliare/Npop;
rawH = dataset.guariti/Npop;
rawE = dataset.deceduti/Npop;

rawT = rawT./max(rawT);
rawD = rawD./max(rawD);
rawH = rawH./max(rawH);
rawE = rawE./max(rawE);

% Matrices for the EKF

std_vec = [5*std(I_data-rawD'), 5*std(I_data-rawD'), std(D_data-rawD'), std(T_data-rawT'), std(H_data-rawH')];
static_pars = opti_lambda.lambda(1);

f_ODE = ParamEvalFunc(f_ODE, static_pars);

ekf = extendedKalmanFilter(@(x) state_transition(f_ODE,x), ...
                           @(x) x(1:5), ...
                           x0);

ekf.StateCovariance = diag([std(S_data), std(I_data), std(D_data), std(T_data), std(H_data), 0.1, 0.1, 0.1, 0.1, 0.1]);
ekf.ProcessNoise = diag([0, 0, 0, 0, 0, 1e-4, 1e-4, 1e-4, 1e-4, 1e-4]);
ekf.MeasurementNoise = diag(std_vec);

x = [x0'; zeros(length(t)-1,length(x0))];
x_lb = [x0'; zeros(length(t)-1,length(x0))];
x_ub = [x0'; zeros(length(t)-1,length(x0))];
Zmeas = [S_data; I_data; D_data; T_data; H_data];

for i = 2:length(t)
    predict(ekf);
    correct(ekf,Zmeas(:,i));
    x(i,:) = ekf.State';
    x_ub(i,:) = ekf.State' + (ones(1,length(x0))*ekf.StateCovariance);
    x_lb(i,:) = ekf.State' - (ones(1,length(x0))*ekf.StateCovariance);
end

%% Plots of the Compartments

close all
set(0,'DefaultFigureWindowStyle','docked');
% Susceptible - S
figure(1)
plot(date, S_data, LineWidth=1.5)
hold on
plot(date, x(:,1), LineWidth=1.5)
hold on
patch([date flip(date)], [x_ub(:,1); flip(x_lb(:,1))], [0.8500, 0.3250, 0.0980], 'FaceAlpha', 0.5, 'EdgeColor', 'none')
ylabel('$\char"0023$ of cases','Interpreter','latex')
title('\textbf{Curve Fitting - Susceptible}','Interpreter','latex')
grid on
legend('Real Data', 'EKF Estimated Data','Interpreter','latex', 'Location','northeast')
xlim([date(1), date(end)])
set(gca, 'TickLabelInterpreter', 'Latex')

% Infected - I
figure(2)
plot(date, I_data, LineWidth=1.5)
hold on
plot(date, x(:,2), LineWidth=1.5)
hold on
patch([date flip(date)], [x_ub(:,2); flip(x_lb(:,2))], [0.8500, 0.3250, 0.0980], 'FaceAlpha', 0.5, 'EdgeColor', 'none')
ylabel('$\char"0023$ of cases','Interpreter','latex')
title('\textbf{Curve Fitting - Infected}','Interpreter','latex')
grid on
legend('Real Data', 'EKF Estimated Data','Interpreter','latex', 'Location','northeast')
xlim([date(1), date(end)])
set(gca, 'TickLabelInterpreter', 'Latex')


% Diagnosed - D
figure(3)
plot(date, D_data, LineWidth=1.5)
hold on
plot(date, x(:,3), LineWidth=1.5)
hold on
patch([date flip(date)], [x_ub(:,3); flip(x_lb(:,3))], [0.8500, 0.3250, 0.0980], 'FaceAlpha', 0.5, 'EdgeColor', 'none')
ylabel('$\char"0023$ of cases','Interpreter','latex')
title('\textbf{Curve Fitting - Diagnosed}','Interpreter','latex')
grid on
legend('Real Data', 'EKF Estimated Data','Interpreter','latex', 'Location','northeast')
xlim([date(1), date(end)])
set(gca, 'TickLabelInterpreter', 'Latex')


% Hospitalised + ICUs - T
figure(4)
plot(date, T_data, LineWidth=1.5)
hold on
plot(date, x(:,4), LineWidth=1.5)
hold on
patch([date flip(date)], [x_ub(:,4); flip(x_lb(:,4))], [0.8500, 0.3250, 0.0980], 'FaceAlpha', 0.5, 'EdgeColor', 'none')
ylabel('$\char"0023$ of cases','Interpreter','latex')
title('\textbf{Curve Fitting - Hospitalised + ICUs }','Interpreter','latex')
grid on
legend('Real Data', 'EKF Estimated Data','Interpreter','latex', 'Location','northeast')
xlim([date(1), date(end)])
set(gca, 'TickLabelInterpreter', 'Latex')


% Healed - H
figure(5)
plot(date, H_data, LineWidth=1.5)
hold on
plot(date, x(:,5), LineWidth=1.5)
hold on
patch([date flip(date)], [x_ub(:,5); flip(x_lb(:,5))], [0.8500, 0.3250, 0.0980], 'FaceAlpha', 0.5, 'EdgeColor', 'none')
ylabel('$\char"0023$ of cases','Interpreter','latex')
title('\textbf{Curve Fitting - Healed}','Interpreter','latex')
grid on
legend('Real Data', 'EKF Estimated Data','Interpreter','latex', 'Location','northeast')
xlim([date(1), date(end)])
set(gca, 'TickLabelInterpreter', 'Latex')


% Coefficient alpha
figure(6)
plot(date, x(:,6), LineWidth=1.5)
hold on
patch([date flip(date)], [x_ub(:,6); flip(x_lb(:,6))], [0.8500, 0.3250, 0.0980], 'FaceAlpha', 0.5, 'EdgeColor', 'none')
ylabel('Coefficients Values','Interpreter','latex')
title('\textbf{Coefficient $\alpha$}','Interpreter','latex')
grid on
legend('$\alpha$','Interpreter','latex', 'Location','northeast')
xlim([date(1), date(end)])
set(gca, 'TickLabelInterpreter', 'Latex')

% Coefficient gamma
figure(7)
plot(date, x(:,7), LineWidth=1.5)
hold on
patch([date flip(date)], [x_ub(:,7); flip(x_lb(:,7))], [0.8500, 0.3250, 0.0980], 'FaceAlpha', 0.5, 'EdgeColor', 'none')
ylabel('Coefficients Values','Interpreter','latex')
title('\textbf{Coefficient $\gamma$}','Interpreter','latex')
grid on
legend('$\gamma$','Interpreter','latex', 'Location','northeast')
xlim([date(1), date(end)])
set(gca, 'TickLabelInterpreter', 'Latex')

% Coefficient delta
figure(8)
plot(date, x(:,8), LineWidth=1.5)
hold on
patch([date flip(date)], [x_ub(:,8); flip(x_lb(:,8))], [0.8500, 0.3250, 0.0980], 'FaceAlpha', 0.5, 'EdgeColor', 'none')
ylabel('Coefficients Values','Interpreter','latex')
title('\textbf{Coefficient $\delta$}','Interpreter','latex')
grid on
legend('$\delta$','Interpreter','latex', 'Location','northeast')
xlim([date(1), date(end)])
set(gca, 'TickLabelInterpreter', 'Latex')

% Coefficient sigma
figure(9)
plot(date, x(:,9), LineWidth=1.5)
hold on
patch([date flip(date)], [x_ub(:,9); flip(x_lb(:,9))], [0.8500, 0.3250, 0.0980], 'FaceAlpha', 0.5, 'EdgeColor', 'none')
ylabel('Coefficients Values','Interpreter','latex')
title('\textbf{Coefficient $\sigma$}','Interpreter','latex')
grid on
legend('$\sigma$','Interpreter','latex', 'Location','northeast')
xlim([date(1), date(end)])
set(gca, 'TickLabelInterpreter', 'Latex')

% Coefficient tau
figure(10)
plot(date, x(:,10), LineWidth=1.5)
hold on
patch([date flip(date)], [x_ub(:,10); flip(x_lb(:,10))], [0.8500, 0.3250, 0.0980], 'FaceAlpha', 0.5, 'EdgeColor', 'none')
ylabel('Coefficients Values','Interpreter','latex')
title('\textbf{Coefficient $\sigma$}','Interpreter','latex')
grid on
legend('$\sigma$','Interpreter','latex', 'Location','northeast')
xlim([date(1), date(end)])
set(gca, 'TickLabelInterpreter', 'Latex')
