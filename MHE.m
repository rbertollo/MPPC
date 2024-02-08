%% Moving Horizon Estimation Trial

clc
clear all
close all
addpath('/Users/marcodelloro/Downloads/casadi-3.6.3-osx64-matlab2018b')
import casadi.*;

load('/Users/marcodelloro/Desktop/Thesis/MSc-Thesis-TUe/Data_Collection/Italian Dataset/SIDTTHE_data_DEF.mat');
set(0,'DefaultFigureWindowStyle','docked');

dataset = readtable('Italy_complete_dataset.xlsx');
dataset = dataset(190:588,:);
dataset.data = datetime(dataset.data, 'InputFormat', 'dd-MMM-yyyy HH:mm:ss');

%%  Data Loading and Initialization

Npop = 59240329; % Total Population of Italy
N = 399;       % MPC horizon length
N_mhe = 21;    % Estimation horizon (3 weeks)
date = SIDTTHE_data{1,1}.date;

I_data = SIDTTHE_data{1,1}.data / Npop;     
D_data = SIDTTHE_data{2,1}.data / Npop;     
T_data = (SIDTTHE_data{3,1}.data + SIDTTHE_data{4,1}.data) / Npop;
H_data = SIDTTHE_data{7,1}.data / Npop;
E_data = SIDTTHE_data{5,1}.data  / Npop;
S_data = ones(length(I_data),1)' - (I_data + D_data + T_data + H_data + E_data );

rawT = dataset.ricoverati + dataset.terapia_intensiva/Npop;
rawD = dataset.isolamento_domiciliare/Npop;
rawH = dataset.guariti/Npop;
rawE = dataset.deceduti/Npop;

rawT = rawT./max(rawT);
rawD = rawD./max(rawD);
rawH = rawH./max(rawH);
rawE = rawE./max(rawE);

y_meas = [S_data; I_data; D_data; T_data; H_data]; % creation of the measurement vector
n_meas = size(y_meas,1);

std_vec = [5*std(I_data-rawD'), 5*std(I_data-rawD'),...
           std(D_data-rawD'), std(T_data-rawT'), std(H_data-rawH')];  % Std measure vector

% "C" measurements matrix from EKF
C = [eye(5), zeros(5,5)];

%% CasADi states initialization

s = casadi.SX.sym('s',1,1); % susceptible population
i = casadi.SX.sym('i',1,1); % infected population
d = casadi.SX.sym('d',1,1); % diagnosed population
t = casadi.SX.sym('t',1,1); % threatned population
h = casadi.SX.sym('h',1,1); % healed population
alp = casadi.SX.sym('alp',1,1);
gma = casadi.SX.sym('gma',1,1);
dlt = casadi.SX.sym('dlt',1,1);
sgm = casadi.SX.sym('sgm',1,1);
tau = casadi.SX.sym('ta',1,1);
v1 = casadi.SX.sym('v1',1,1);
v2 = casadi.SX.sym('v2',1,1);
v3 = casadi.SX.sym('v3',1,1);
v4 = casadi.SX.sym('v4',1,1);
v5 = casadi.SX.sym('v5',1,1);

lambda = 0.1;

x = [s; i; d; t; h; alp; gma; dlt; sgm; tau];


n_states = length(x);

%% Do the state update by hand and close the loop (multiple shooting)

opti = casadi.Opti();
X = opti.variable(n_states, N_mhe); % this represents all the states subjected to an actual dynamics
p = opti.parameter(5,1);
v = opti.variable(5,N_mhe-1);

eqns = [    -x(1) * (x(6) * x(2)) + lambda * x(2);...
             x(1) * (x(6) * x(2)) - (x(7)+lambda) * x(2);...
             x(2) * x(7) - x(3) * (lambda + x(8));...
             x(8) * x(3) - (x(10) + x(9))*x(4);...
             lambda * x(3) + x(4) * x(9);...
             0;...
             0;...
             0;...
             0;...
             0          ];

f = casadi.Function('f',{x},{eqns}); % Model Dynamics

% Computation of the dynamics constraint
Ts = 1;

for k=1:N_mhe-1
   
    % Runge-Kutta 4 integration
    k1 = f(X(:,k));
    k2 = f(X(:,k)+Ts/2*k1);
    k3 = f(X(:,k)+Ts/2*k2);
    k4 = f(X(:,k)+Ts*k3);
    x_plus = Ts/6*(k1+2*k2+2*k3+k4);
    opti.subject_to(X(:,k+1)==X(:,k) + x_plus + [v(:,k); zeros(5,1)]); % close the gaps - dynamics consstraint

end

%% Optimization algorithm

matrix_par = [];
matrix_sts = [];
stdevs = std(y_meas,1,2);

opti.subject_to(X(:)>0);
c = 0.005;

opti.subject_to(-stdevs(1)*c < v(1,:) < stdevs(1)*c);
opti.subject_to(-stdevs(2)*c < v(2,:) < stdevs(2)*c);
opti.subject_to(-stdevs(3)*c < v(3,:) < stdevs(3)*c);
opti.subject_to(-stdevs(4)*c < v(4,:) < stdevs(4)*c);
opti.subject_to(-stdevs(5)*c < v(5,:) < stdevs(5)*c);

for ii=1:N_mhe
    opti.subject_to( 1 - sum(X(1:5,ii)) >= 0  );
end

opti.subject_to( 0.05 <= X(6,:)<= 0.8 ); % alpha constr.
opti.subject_to( 0.005 <= X(7,:) <= 0.6 ); % gamma constr.
opti.subject_to( 1e-4 <= X(8,:) <= 0.6 ); % delta constr.
opti.subject_to( 1e-4 <= X(9,:) <= 1 ); % sigma constr.
opti.subject_to( 1e-4 <= X(10,:) <= 0.5 ); % tau constr.

% options on solver
opts.print_time = 0;
opts.ipopt.print_level = 0;
opti.solver('ipopt', opts);

xtilde = y_meas(1:5,1);    
ptilde = [ 0.4; 0.3; 0.05; 0.1; 0.35 ];

for k = N_mhe:1:N-N_mhe

    Q = diag([ 0.5 0.5 5 5 5]);
    obj = []; % initialization of the object for every iteration 
    obj = sumsqr(X(1:5,:) - xtilde) + 10*sumsqr(X(6:10,:) - ptilde) + sumsqr( Q*(y_meas(:,k-N_mhe+1:k) - X(1:5,:))./std(y_meas,1,2) ) + sumsqr(v);   
    opti.minimize(obj);

    sol = opti.solve();
    
    % Save of the optimization values
    params_dyn = opti.value(X(6:end,end))';
    states_dyn = opti.value(X(1:5,end))';

    % Append the current row to the main table
    matrix_par = [matrix_par; params_dyn];
    matrix_sts = [matrix_sts; states_dyn];
    ptilde = opti.value(X(6:10,1));
    xtilde = opti.value(X(1:5,2));

end

%% Plot of the results

% Results saving in a results struct
column_names_par = {'alpha', 'gamma', 'delta','sigma', 'tau'};
column_names_sts = {'S', 'I', 'D', 'T', 'H'};

table_par = array2table(matrix_par, 'VariableNames', column_names_par);
table_sts = array2table(matrix_sts, 'VariableNames', column_names_sts);

results.par = table_par;
results.sts = table_sts;

%% Additional Intersting plots,

policy_idx = [1 40 69 116 141 190 242 294 399];   % Values taken from italian policies applied for the pandemic

for ii=1:length(policy_idx)
    policy_dates(ii) = [date(policy_idx(ii))];
end

customColors2 = {   [1, 0.6, 0.2],
                    [1, 0.3, 0.05],
                    [1, 0.2, 0.05],
                    [0.8, 0.2, 0.1],
                    [1, 0.2, 0.05],            
                    [0.8, 0.2, 0.1],
                    [1, 0.3, 0.05],
                    [1, 0.6, 0.2]
                };

for ii = 1:length(policy_dates)-1

    area.x(ii, :) = [policy_dates(ii) policy_dates(ii) policy_dates(ii+1) policy_dates(ii+1)];
    area.y_alpha(ii, :) = [min(results.par.alpha)*0.5 max(results.par.alpha)*1.05 max(results.par.alpha)*1.05 min(results.par.alpha)*0.5];
end

% States
figure(1)
plt=plot(date(N_mhe:N-N_mhe),results.sts.S, 's-', 'LineWidth', 1, 'MarkerSize', 5);
plt.MarkerFaceColor = plt.Color;
hold on
plot(date(N_mhe:N-N_mhe), S_data(N_mhe:N-N_mhe), LineWidth=1.5)
xlim([date(1+N_mhe), date(end-N_mhe)])
title('\textbf{S}','Interpreter','latex')

figure(2)
plt=plot(date(N_mhe:N-N_mhe),results.sts.I, 's-', 'LineWidth', 1, 'MarkerSize', 5);
plt.MarkerFaceColor = plt.Color;
hold on 
plot(date(N_mhe:N-N_mhe), I_data(N_mhe:N-N_mhe), LineWidth=1.5)
xlim([date(1+N_mhe), date(end-N_mhe)])
title('\textbf{I}','Interpreter','latex')

figure(3)
plt=plot(date(N_mhe:N-N_mhe),results.sts.D, 's-', 'LineWidth', 1, 'MarkerSize', 5);
plt.MarkerFaceColor = plt.Color;
hold on 
plot(date(N_mhe:N-N_mhe), D_data(N_mhe:N-N_mhe), LineWidth=1.5)
xlim([date(1+N_mhe), date(end-N_mhe)])
title('\textbf{D}','Interpreter','latex')

figure(4)
plt=plot(date(N_mhe:N-N_mhe),results.sts.T, 's-', 'LineWidth', 1, 'MarkerSize', 5);
plt.MarkerFaceColor = plt.Color;
hold on 
plot(date(N_mhe:N-N_mhe), T_data(N_mhe:N-N_mhe), LineWidth=1.5)
xlim([date(1+N_mhe), date(end-N_mhe)])
title('\textbf{T}','Interpreter','latex')

figure(5)
plt=plot(date(N_mhe:N-N_mhe),results.sts.H, 's-', 'LineWidth', 1, 'MarkerSize', 5);
plt.MarkerFaceColor = plt.Color;
hold on 
plot(date(N_mhe:N-N_mhe), H_data(N_mhe:N-N_mhe), LineWidth=1.5)
xlim([date(1+N_mhe), date(end-N_mhe)])
title('\textbf{H}','Interpreter','latex')

% Figure of the alpha trend related to policy In italy
figure(10)
for ii = 1:length(policy_dates)-1
    fill(area.x(ii, :) ,area.y_alpha(1, :), customColors2{ii,1} ,'FaceAlpha',.5,'EdgeColor', 'none','HandleVisibility', 'off')
    hold on
    xline(policy_dates(ii),":",'HandleVisibility', 'off')
    hold on 
end
hold on
plot(date(N_mhe:N-N_mhe), results.par.alpha,'k','LineWidth',1.5, 'DisplayName', '$\alpha$')
ylabel('Coefficients Values','Interpreter','latex')
title('\textbf{$\alpha$ coefficient}','Interpreter','latex')
grid on
legend('Interpreter','latex','location','southeast')
xlim([date(1+N_mhe), date(end-N_mhe)])
ylim([min(results.par.alpha)*0.5, max(results.par.alpha)*1.05])
set(gca, 'TickLabelInterpreter', 'Latex')

figure(6)
plot(date(N_mhe:N-N_mhe),results.par.gamma, 'LineWidth', 1.5, 'Color','k')
xlim([date(1+N_mhe), date(end-N_mhe)])
title('\textbf{$\gamma$}','Interpreter','latex')

figure(7)
plot(date(N_mhe:N-N_mhe),results.par.delta, 'LineWidth', 1.5, 'Color','k')
xlim([date(1+N_mhe), date(end-N_mhe)])
title('\textbf{$\delta$}','Interpreter','latex')

figure(8)
plot(date(N_mhe:N-N_mhe),results.par.sigma, 'LineWidth', 1.5, 'Color','k')
xlim([date(1+N_mhe), date(end-N_mhe)])
title('\textbf{$\sigma$}','Interpreter','latex')

figure(9)
plot(date(N_mhe:N-N_mhe),results.par.tau, 'LineWidth', 1.5, 'Color','k')
xlim([date(1+N_mhe), date(end-N_mhe)])
title('\textbf{$\tau$}','Interpreter','latex')