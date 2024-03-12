%% Moving Horizon Estimation Trial

clc
clear all
close all
addpath('/Users/marcodelloro/Downloads/casadi-3.6.3-osx64-matlab2018b')
import casadi.*;

load('/Users/marcodelloro/Desktop/Thesis/MSc-Thesis-TUe/Data_Collection/Italian Dataset/SIDTTHE_data_DEF.mat');
% set(0,'DefaultFigureWindowStyle','docked');

dataset = readtable('Italy_complete_dataset.xlsx');
dataset = dataset(190:588,:);
dataset.data = datetime(dataset.data, 'InputFormat', 'dd-MMM-yyyy HH:mm:ss');
load('/Users/marcodelloro/Desktop/Thesis/MSc-Thesis-TUe/Data_Collection/Italian Dataset/Var_Infected.mat')

%%  Data Loading and Initialization

Npop = 59240329; % Total Population of Italy
N = 399;       % MPC horizon length
N_mhe = 21;    % Estimation horizon (3 weeks)
date = SIDTTHE_data{1,1}.date;

I_data = SIDTTHE_data{1,1}.data / Npop;     
D_data = SIDTTHE_data{2,1}.data / Npop;     
T_data = (SIDTTHE_data{3,1}.data + SIDTTHE_data{4,1}.data) / Npop;
H_data = SIDTTHE_data{7,1}.data / Npop;
H_dataAug = SIDTTHE_data{6,1}.data / Npop;
E_data = SIDTTHE_data{5,1}.data  / Npop;
S_data = ones(length(I_data),1)' - (I_data + D_data + T_data + H_dataAug + E_data );

rawT = dataset.ricoverati + dataset.terapia_intensiva/Npop;
rawD = dataset.isolamento_domiciliare/Npop;
rawH = dataset.guariti/Npop;
rawE = dataset.deceduti/Npop;

rawT = rawT./max(rawT);
rawD = rawD./max(rawD);
rawH = rawH./max(rawH);
rawE = rawE./max(rawE);

y_meas = [S_data; I_data; D_data; T_data; H_dataAug; E_data]; % creation of the measurement vector
n_meas = size(y_meas,1);

std_vec = [5*std(I_data-rawD'), 5*std(I_data-rawD'),...
           std(D_data-rawD'), std(T_data-rawT'), std(H_data-rawH'), std(E_data-rawE')];  % Std measure vector

% "C" measurements matrix from EKF
C = [eye(6), zeros(6,5)];

%% CasADi states initialization

s = casadi.SX.sym('s',1,1); % susceptible population
i = casadi.SX.sym('i',1,1); % infected population
d = casadi.SX.sym('d',1,1); % diagnosed population
t = casadi.SX.sym('t',1,1); % threatned population
h = casadi.SX.sym('h',1,1); % healed population
e = casadi.SX.sym('e',1,1); % expired population
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

x = [s; i; d; t; h; e; alp; gma; dlt; sgm; tau];

n_states = length(x);

%% Do the state update by hand and close the loop (multiple shooting)

opti = casadi.Opti();
X = opti.variable(n_states, N_mhe); % this represents all the states subjected to an actual dynamics
p = opti.parameter(5,1);
v = opti.variable(6,N_mhe-1);

eqns = [    -x(1) * (x(7) * x(2));...
             x(1) * (x(7) * x(2)) - (x(8)+lambda) * x(2);...
             x(2) * x(8) - x(3) * (lambda + x(9));...
             x(9) * x(3) - ( x(11) + x(10) )*x(4);...
             lambda * x(3) + x(4) * x(10) + lambda * x(2);...
             x(11) * x(4);...
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
    opti.subject_to(X(:,k+1)==X(:,k) + x_plus + [v(:,k); zeros(5,1)]); % close the gaps - dynamics constraint
    
end

%% Optimization algorithm

matrix_par = [];
matrix_sts = [];
matrix_szero = [];

stdevs = std(y_meas,1,2);

opti.subject_to(X(:)>0);
c = 0.005;

opti.subject_to(-stdevs(1)*c < v(1,:) < stdevs(1)*c);
opti.subject_to(-stdevs(2)*c < v(2,:) < stdevs(2)*c);
opti.subject_to(-stdevs(3)*c < v(3,:) < stdevs(3)*c);
opti.subject_to(-stdevs(4)*c < v(4,:) < stdevs(4)*c);
opti.subject_to(-stdevs(5)*c < v(5,:) < stdevs(5)*c);
opti.subject_to(-stdevs(6)*c < v(6,:) < stdevs(6)*c);

for ii=1:N_mhe
    % opti.subject_to( 1 - sum(X(1:6,ii)) >= 0  );
    opti.subject_to( sum(X(1:6,ii)) == 1 )
end

opti.subject_to( 0.05 <= X(7,:)<= 0.8 ); % alpha constr.
opti.subject_to( 0.005 <= X(8,:) <= 0.6 ); % gamma constr.
opti.subject_to( 1e-4 <= X(9,:) <= 0.6 ); % delta constr.
opti.subject_to( 1e-4 <= X(10,:) <= 0.5 ); % sigma constr.
opti.subject_to( 1e-4 <= X(11,:) <= 0.5 ); % tau constr.

% options on solver
opts.print_time = 0;
opts.ipopt.print_level = 0;
opti.solver('ipopt', opts);

xtilde = y_meas(1:6,1);    
ptilde = [0.25; 0.12; 0.01; 0.02; 0.02];
dyn = [xtilde' ptilde' ];
Pi_min = diag([std(S_data), std(I_data), std(D_data), std(T_data), std(H_data), std(E_data) 0.01, 0.01, 0.01, 0.01, 0.01]);

% cost 1 = X(1:5,:) - xtilde --> minimize difference between previous
%                                estimation and current estimation (states)
% cost 2 = X(6:10,:) - ptilde --> minimize difference between previous
%                                estimation and current estimation (parameters)      
% cost 3 = y_meas - X(1:5,:) - ptilde --> minimize difference measure and
%                                         estimate, at curent state
% cost 4 = xdot - f(x) == v(t) --> minimize process noise at the current time instant

% Names of the states
column_names_sts = {'S', 'I', 'D', 'T', 'H', 'E'};

%Initialization of ful state matri
iterationIndex = 1;

for k = N_mhe:1:N-N_mhe
    % matrices for EKF cov. update
    [A,C,G] = getMatricesEKF(dyn);
    Q = diag([1e-4, 1e-4, 1e-4, 1e-4, 1e-4, 1e-4, 0.1, 0.1, 0.1, 0.1, 0.1]);
    R = [0.0387	0.0387 0.00570 0.0080 0.00076 0.006];
    Z1 = diag([1 1 5 5 2 1]);
    Z2 = diag([1.5 25 40 15 10]);
    Pi = G*Q*G' + A*Pi_min*A' - A*Pi_min*C'*inv(R + C*Pi_min*C')*C*Pi_min*A';
    Zinitial_cond = diag([1 1 1 1 1 1]);
    cost1 = [X(1:6,1) - xtilde; X(7:11,1) - ptilde]'*inv(Pi)*[X(1:6,1) - xtilde; X(7:11,1) - ptilde];
 
    obj = []; % initialization of the object for every iteration 
    obj = sumsqr(Zinitial_cond*(X(1:6,:) - xtilde)) + 100*sumsqr(Z2*(X(7:11,:) - ptilde)) + sumsqr( Z1*(y_meas(:,k-N_mhe+1:k) - X(1:6,:))./std(y_meas,1,2) ) + 100*sumsqr(v); 
    % obj = 100*cost1 + sumsqr( Z1*(y_meas(:,k-N_mhe+1:k) - X(1:6,:))./std(y_meas,1,2) ) + 100*sumsqr(v);   
    opti.minimize(obj);

    sol = opti.solve();
    dyn = opti.value(X(:,end));
    
    % Save of the optimization values
    params_dyn = opti.value(X(7:end,end))';
    states_dyn = opti.value(X(1:6,end))';
    states_zero = opti.value(X(1:6,1))';
    
    % CURRENT STATE SAVE
    currentState = opti.value(X(1:6,:))';
    currentTable = array2table(currentState, 'VariableNames', column_names_sts);
    stateTables{iterationIndex} = currentTable; % Saving all the tables in the cell struct
    iterationIndex = iterationIndex + 1;

    % Append the current row to the main table
    matrix_par = [matrix_par; params_dyn];
    matrix_sts = [matrix_sts; states_dyn];
    matrix_szero = [matrix_szero; states_zero];
   
    ptilde = opti.value(X(7:11,1));
    xtilde = opti.value(X(1:6,2));
end

%% Plot of the results

% Results saving in a results struct
column_names_par = {'alpha', 'gamma', 'delta','sigma', 'tau'};

table_par = array2table(matrix_par, 'VariableNames', column_names_par);
table_sts = array2table(matrix_sts, 'VariableNames', column_names_sts);
table_sts_zero = array2table(matrix_szero, 'VariableNames', column_names_sts);

results.par = table_par;
results.sts = table_sts;
results.stszero = table_sts_zero;
results.FullStates = stateTables;

% save the results in a CSV to import in into Python
path = '/Users/marcodelloro/Desktop/Thesis/MSc-Thesis-TUe/Forecasting/TrialPython';
filePathPar = fullfile(path, 'table_par2.csv');
writetable(results.par, filePathPar);

save('results2.mat', 'results');

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
    area.y_alpha(ii, :) = [0 max(results.par.alpha)*1.5 max(results.par.alpha)*1.5 0];
end

% States
figure(1)
scatter(date(N_mhe:N-N_mhe), S_data(N_mhe:N-N_mhe),20,'filled')
hold on 
plot(date(N_mhe:N-N_mhe),results.sts.S, 'LineWidth', 1.5, 'MarkerSize', 5);
xlim([date(1+N_mhe), date(end-N_mhe)])
title('\textbf{S}','Interpreter','latex')
grid on
set(gca, 'TickLabelInterpreter', 'Latex')
lgd = legend('Real Data', 'MHE Fitted Data','Interpreter','latex','location','northeast');
lgd.FontSize = 18;

figure(2)
sct = scatter(date(N_mhe:N-N_mhe), I_data(N_mhe:N-N_mhe),100,'filled','MarkerEdgeAlpha', 0.3,'MarkerFaceAlpha', 0.3); % Set face transparency
colorsct=[0, 0.4470, 0.7410];
hold on 
fill(var_area.x,var_area.Ivar/Npop, colorsct, 'FaceAlpha', .1, 'EdgeColor', 'none');
hold on
plot(date(N_mhe:N-N_mhe),results.sts.I, 'LineWidth', 2, 'MarkerSize', 5,'Color',[0.8500 0.3250 0.0980]);
xlim([date(1+N_mhe), date(end-N_mhe)])
% title('\textbf{\textit{I}-Infected individuals}','Interpreter','latex')
yax = ylabel('Normalized Population','Interpreter','latex');
% yax.FontSize = 14;
grid on
set(gca, 'TickLabelInterpreter', 'Latex')
lgd = legend('Estimated Data', '95\% confidence interval','MHE Fitted Data','Interpreter','latex','location','northeast');
% lgd.FontSize = 14;
box on 

figure(3)
scatter(date(N_mhe:N-N_mhe), D_data(N_mhe:N-N_mhe),100,'filled','MarkerEdgeAlpha', 0.8,'MarkerFaceAlpha', 0.8)
hold on 
plot(date(N_mhe:N-N_mhe),results.sts.D, 'LineWidth', 2, 'MarkerSize', 5);
xlim([date(1+N_mhe), date(end-N_mhe)])
% title('\textbf{\textit{D} - Detected individuals}','Interpreter','latex')
grid on
set(gca, 'TickLabelInterpreter', 'Latex')
yax = ylabel('Normalized Population','Interpreter','latex');
% yax.FontSize = 14;
lgd = legend('Real Data', 'MHE Fitted Data','Interpreter','latex','location','northeast');
% lgd.FontSize = 14;
box on 

figure(4)
scatter(date(N_mhe:N-N_mhe), T_data(N_mhe:N-N_mhe),100,'filled','MarkerEdgeAlpha', 0.8,'MarkerFaceAlpha', 0.8)
hold on 
plot(date(N_mhe:N-N_mhe),results.sts.T, 'LineWidth', 2, 'MarkerSize', 5);
xlim([date(1+N_mhe), date(end-N_mhe)])
% title('\textbf{\textit{T} - Hospitalised and ICUs individuals}','Interpreter','latex')
grid on
set(gca, 'TickLabelInterpreter', 'Latex')
yax = ylabel('Normalized Population','Interpreter','latex');
box on 
% yax.FontSize = 14;
% lgd = legend('Real Data', 'MHE Fitted Data','Interpreter','latex','location','northeast');
% lgd.FontSize = 18;

figure(5)
scatter(date(N_mhe:N-N_mhe), H_data(N_mhe:N-N_mhe),100,'filled','MarkerEdgeAlpha', 0.8,'MarkerFaceAlpha', 0.8)
hold on 
% scatter(date(N_mhe:N-N_mhe), H_dataAug(N_mhe:N-N_mhe),20,'filled')
% hold on 
plot(date(N_mhe:N-N_mhe),results.sts.H, 'LineWidth', 2, 'MarkerSize', 5);
% hold on 
% plot(date(N_mhe:N-N_mhe),results.stszero.H, 'LineWidth', 1.5, 'MarkerSize', 5);
xlim([date(1+N_mhe), date(end-N_mhe)])
% title('\textbf{\textit{H} - Healed individuals}','Interpreter','latex')
grid on
set(gca, 'TickLabelInterpreter', 'Latex')
yax = ylabel('Normalized Population','Interpreter','latex');
box on 
% yax.FontSize = 14;
% lgd = legend('Real Data', 'Augmented Data','MHE Fitted Data','MHE $x_0$', 'Interpreter','latex','location','northeast');
% lgd.FontSize = 18;

figure(106)
scatter(date(N_mhe:N-N_mhe), E_data(N_mhe:N-N_mhe),100,'filled','MarkerEdgeAlpha', 0.8,'MarkerFaceAlpha', 0.8)
hold on 
plot(date(N_mhe:N-N_mhe),results.sts.E, 'LineWidth', 2, 'MarkerSize', 5);
hold on 
% plot(date(N_mhe:N-N_mhe),results.stszero.E, 'LineWidth', 1.5, 'MarkerSize', 5);
xlim([date(1+N_mhe), date(end-N_mhe)])
% title('\textbf{\textit{D} - Expired individuals}','Interpreter','latex')
grid on
set(gca, 'TickLabelInterpreter', 'Latex')
yax = ylabel('Normalized Population','Interpreter','latex');
box on 
% yax.FontSize = 14;
% lgd = legend('Real Data', 'MHE Fitted Data','Fitted $x_0$','Interpreter','latex','location','northeast');
% lgd.FontSize = 18;


% Figure of the alpha trend related to policy In italy
figure(10)
for ii = 1:length(policy_dates)-1
    handleVisibilityValue = 'off';
    if ii >= 1 && ii <= 4
        handleVisibilityValue = 'on';
    end
    
    fill(area.x(ii, :) ,area.y_alpha(1, :), customColors2{ii,1} ,...
        'FaceAlpha',.5,'EdgeColor', 'none','HandleVisibility', handleVisibilityValue)
    hold on
    xline(policy_dates(ii),":",'HandleVisibility', 'off')
    hold on 
end
hold on
plot(date(N_mhe:N-N_mhe), results.par.alpha,'k','LineWidth',2, 'DisplayName', '$\alpha$','HandleVisibility', 'off')
yax = ylabel('Coefficient Value','Interpreter','latex');
% yax.FontSize = 14;
% title('\textbf{$\alpha$ coefficient}','Interpreter','latex')
grid on
lgd = legend('Mild Restrictions', 'NPIs \& Social Limitations','Curfew','Total Lockdown','Interpreter','latex','location','northeast');
% lgd.FontSize = 18;
title(lgd, '\textbf{Policy Level}');
xlim([date(1+N_mhe), date(end-N_mhe)])
ylim([0, max(results.par.alpha)*1.5])
set(gca, 'TickLabelInterpreter', 'Latex')


figure(6)
plot(date(N_mhe:N-N_mhe),results.par.gamma, 'LineWidth', 2, 'Color','k')
xlim([date(1+N_mhe), date(end-N_mhe)])
ylim([0, max(results.par.gamma)*1.5])
grid on
set(gca, 'TickLabelInterpreter', 'Latex')
yax = ylabel('Coefficient Value','Interpreter','latex');
% yax.FontSize = 14;
% title('\textbf{$\gamma$}','Interpreter','latex')

figure(7)
plot(date(N_mhe:N-N_mhe),results.par.delta, 'LineWidth', 2, 'Color','k')
xlim([date(1+N_mhe), date(end-N_mhe)])
ylim([0, max(results.par.delta)*1.2])
grid on
set(gca, 'TickLabelInterpreter', 'Latex')
yax = ylabel('Coefficient Value','Interpreter','latex');
% yax.FontSize = 14;
% title('\textbf{$\delta$}','Interpreter','latex')

figure(8)
plot(date(N_mhe:N-N_mhe),results.par.sigma, 'LineWidth', 2, 'Color','k')
xlim([date(1+N_mhe), date(end-N_mhe)])
ylim([0, max(results.par.sigma)*1.5])
grid on
set(gca, 'TickLabelInterpreter', 'Latex')
% title('\textbf{$\sigma$}','Interpreter','latex')
yax = ylabel('Coefficient Value','Interpreter','latex');
% yax.FontSize = 14;

figure(9)
plot(date(N_mhe:N-N_mhe),results.par.tau, 'LineWidth', 2, 'Color','k')
xlim([date(1+N_mhe), date(end-N_mhe)])
ylim([0, max(results.par.tau)*1.5])
grid on
set(gca, 'TickLabelInterpreter', 'Latex')
yax = ylabel('Coefficient Value','Interpreter','latex');
% yax.FontSize = 14;
% title('\textbf{$\tau$}','Interpreter','latex')

mape_error = struct();
struct_idx = {'S', 'I', 'D', 'T', 'H', 'E'};
y_meas2 = [S_data; I_data; D_data; T_data; H_data; E_data];

% Iterate over y_meas
for i = 1:size(y_meas2,1)
    y_true = y_meas2(i,:);
    y_true = y_true(N_mhe:N-N_mhe);
    y_pred = results.sts.(struct_idx{i})';
    mape = mapeFunc(y_true, y_pred);        % MAPE formula

    mape_error.(struct_idx{i}) = mape;
    rmse = sqrt(mean((y_true - y_pred).^2));  % RMSE formula
    rmse_error.(struct_idx{i}) = rmse;  % Store the RMSE in the structure
end
%%
N_mhe = 21;    % Estimation horizon (3 weeks)
StartHorizon = 10*N_mhe; % we are here

mape_errorFore = struct();
struct_idx = {'S', 'I', 'D', 'T', 'H', 'E'};
y_meas2 = [S_data; I_data; D_data; T_data; H_data; E_data];

% Iterate over y_meas
weekAfter = [0,7,14,21,28,35,42,49];

for jj = 1:8
    StartHorizon = 10*N_mhe + weekAfter(jj);
    for i = 1:size(y_meas2,1)
        y_true = y_meas2(i,:);
        y_true = y_true(N_mhe:StartHorizon+N_mhe-1);  % Adjusted indexing to match MATLAB's 1-based indexing
        y_pred = results.sts.(struct_idx{i})(1:StartHorizon)';
        mape = mapeFunc(y_true, y_pred);  % MAPE formula

        if ~isfield(mape_errorFore, ['weekAfter' num2str(weekAfter(jj))])
            mape_errorFore.(['weekAfter' num2str(weekAfter(jj))]) = struct();
        end
        % Store the mape for the current compartment and weekAfter index
        mape_errorFore.(['weekAfter' num2str(weekAfter(jj))]).(struct_idx{i}) = mape;
    end
end