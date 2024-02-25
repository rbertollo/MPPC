%% Moving Horizon Estimation + Forecast
% Data Loading

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

load("results.mat")

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

y_meas = [S_data; I_data; D_data; T_data; H_dataAug; E_data]; % creation of the measurement vector

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

% Equation for the ODEs simulation
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

%% Simulation in the future of the pandemic spread using the model - CONSTANT PARAMETERS

stepsaheadVec = [21, 35];
N_mhe = 21;
StartHorizon = 4*N_mhe; % ATTENTION !! If this is changed, must be changed also in Python script
weeksIn = [7 14 21 28 35];
weeksName = [0 1 2 3 4 5]; % Vector to save better the results 
FitHorizon = [  StartHorizon; StartHorizon+weeksIn(1); StartHorizon+weeksIn(2);...
                StartHorizon+weeksIn(3); StartHorizon+weeksIn(4); StartHorizon+weeksIn(5) ];

for jj = 1:length(FitHorizon)
    for i = 1:length(stepsaheadVec)
        stepsahead = stepsaheadVec(i);
        
        X0 = [  y_meas(1,FitHorizon(jj));... S starting cond.
                y_meas(2,FitHorizon(jj));... I starting cond.
                y_meas(3,FitHorizon(jj));... D starting cond.
                y_meas(4,FitHorizon(jj));... T starting cond.
                y_meas(5,FitHorizon(jj));... H starting cond.
                y_meas(6,FitHorizon(jj));... E starting cond.
                results.par.alpha(FitHorizon(jj));...
                results.par.gamma(FitHorizon(jj));...
                results.par.delta(FitHorizon(jj));...
                results.par.sigma(FitHorizon(jj));...
                results.par.tau(FitHorizon(jj))          ];
        
        % Initialize the state matrix simulation
        X_sim = zeros(n_states, stepsahead);
        X_sim(:,1) = X0;
        Ts = 1;
    
        for kk = 1:stepsahead-1
            % Runge-Kutta 4 integration
            k1 = full(f(X_sim(:,kk)));
            k2 = full(f(X_sim(:,kk) + Ts/2*k1));
            k3 = full(f(X_sim(:,kk) + Ts/2*k2));
            k4 = full(f(X_sim(:,kk) + Ts*k3));
            x_plus = Ts/6*(k1 + 2*k2 + 2*k3 + k4);
            X_sim(:,kk+1) = X_sim(:,kk) + full(x_plus);
        end
       
        fieldName = sprintf('ForeWk%d_steps%d', weeksName(jj), stepsahead);
        Xsim.(fieldName) = X_sim;
        
    end
end

%% Data Loading & Prep from Python simulation

ForeWk0 = readtable('/Users/marcodelloro/Desktop/Thesis/MSc-Thesis-TUe/Forecasting/TrialPython/Forecasts21d_CI.xlsx', 'Sheet', 'Horizon_0_days');
ForeWk1 = readtable('/Users/marcodelloro/Desktop/Thesis/MSc-Thesis-TUe/Forecasting/TrialPython/Forecasts21d_CI.xlsx', 'Sheet', 'Horizon_7_days');
ForeWk2 = readtable('/Users/marcodelloro/Desktop/Thesis/MSc-Thesis-TUe/Forecasting/TrialPython/Forecasts21d_CI.xlsx', 'Sheet', 'Horizon_14_days');
ForeWk3 = readtable('/Users/marcodelloro/Desktop/Thesis/MSc-Thesis-TUe/Forecasting/TrialPython/Forecasts21d_CI.xlsx', 'Sheet', 'Horizon_21_days');
ForeWk4 = readtable('/Users/marcodelloro/Desktop/Thesis/MSc-Thesis-TUe/Forecasting/TrialPython/Forecasts21d_CI.xlsx', 'Sheet', 'Horizon_28_days');
ForeWk5 = readtable('/Users/marcodelloro/Desktop/Thesis/MSc-Thesis-TUe/Forecasting/TrialPython/Forecasts21d_CI.xlsx', 'Sheet', 'Horizon_35_days');

% Adjusting al the data and reorganize them in Struct - Easier to handle 
columnHolt = {'alpha_yhat','alpha_yhat_lower','alpha_yhat_upper','gamma_yhat','gamma_yhat_lower','gamma_yhat_upper','delta_yhat','delta_yhat_lower','delta_yhat_upper'...
              'sigma_yhat','sigma_yhat_lower','sigma_yhat_upper','tau_yhat','tau_yhat_lower','tau_yhat_upper'};

columnProph = {'alpha_Prophet_yhat','alpha_Prophet_yhat_lower','alpha_Prophet_yhat_upper','gamma_Prophet_yhat','gamma_Prophet_yhat_lower','gamma_Prophet_yhat_upper',...
              'delta_Prophet_yhat','delta_Prophet_yhat_lower','delta_Prophet_yhat_upper','sigma_Prophet_yhat','sigma_Prophet_yhat_lower','sigma_Prophet_yhat_upper',...
              'tau_Prophet_yhat','tau_Prophet_yhat_lower','tau_Prophet_yhat_upper'};

% Struct for Holt forecasting
HoltFore.wk0 = ForeWk0(:, columnHolt);
HoltFore.wk1 = ForeWk1(:, columnHolt);
HoltFore.wk2 = ForeWk2(:, columnHolt);
HoltFore.wk3 = ForeWk3(:, columnHolt);
HoltFore.wk4 = ForeWk4(:, columnHolt);
HoltFore.wk5 = ForeWk5(:, columnHolt);

% Struct for Prophet forecasting

ProphFore.wk0 = ForeWk0(:, columnProph);
ProphFore.wk1 = ForeWk1(:, columnProph);
ProphFore.wk2 = ForeWk2(:, columnProph);
ProphFore.wk3 = ForeWk3(:, columnProph);
ProphFore.wk4 = ForeWk4(:, columnProph);
ProphFore.wk5 = ForeWk5(:, columnProph);

s = casadi.SX.sym('s',1,1); 
i = casadi.SX.sym('i',1,1);
d = casadi.SX.sym('d',1,1); 
t = casadi.SX.sym('t',1,1); 
h = casadi.SX.sym('h',1,1); 
e = casadi.SX.sym('e',1,1); 
m = [s; i; d; t; h; e];
m_states = length(m);


%% Simulation in the future of the pandemic spread using the model - Holt's estimated parameters

X_m = zeros(m_states, stepsaheadVec(1));
paramSets = {'yhat', 'yhat_lower', 'yhat_upper'};

for jj = 1:length(FitHorizon)
    
    % X0 = [  y_meas(1,FitHorizon(jj));... S starting cond.
    %         y_meas(2,FitHorizon(jj));... I starting cond.
    %         y_meas(3,FitHorizon(jj));... D starting cond.
    %         y_meas(4,FitHorizon(jj));... T starting cond.
    %         y_meas(5,FitHorizon(jj));... H starting cond.
    %         y_meas(6,FitHorizon(jj))]; % E starting cond.

    % starting condition MHE
    X0 = [  results.sts.S(FitHorizon(jj));... S starting cond.
            results.sts.I(FitHorizon(jj));... I starting cond.
            results.sts.D(FitHorizon(jj));... D starting cond.
            results.sts.T(FitHorizon(jj));... T starting cond.
            results.sts.H(FitHorizon(jj));... H starting cond.
            results.sts.E(FitHorizon(jj))]; % E starting cond.
    
    fieldName = sprintf('wk%d', jj-1);
    HoltParams = getfield(HoltFore, fieldName);

    for ii = 1:length(paramSets)
        X_m(:,1) = X0;

        for k=1:stepsaheadVec(1)-1
            paramSuffix = paramSets{ii};

            alpha = HoltParams.(['alpha_', paramSuffix])(k);
            gamma = HoltParams.(['gamma_', paramSuffix])(k);
            delta = HoltParams.(['delta_', paramSuffix])(k);
            sigma = HoltParams.(['sigma_', paramSuffix])(k);
            tau   = HoltParams.(['tau_', paramSuffix])(k);

             eqns2 = [   -m(1) * (alpha * m(2));...
                          m(1) * (alpha * m(2)) - (gamma+lambda) * m(2);...
                          m(2) * gamma - m(3) * (lambda + delta);...
                          delta * m(3) - ( tau + sigma )*m(4);...
                          lambda * m(3) + m(4) * sigma + lambda * m(2);...
                          tau * m(4)     ];

        f2 = casadi.Function('f2', {m}, {eqns2});

        k1 = f2(X_m(:,k));
        k2 = f2(X_m(:,k)+Ts/2*k1);
        k3 = f2(X_m(:,k)+Ts/2*k2);
        k4 = f2(X_m(:,k)+Ts*k3);
        Xm_plus = Ts/6*(k1+2*k2+2*k3+k4);
        X_m(:,k+1) = X_m(:,k) + full(Xm_plus);

        end

        HoltWintSim.(fieldName).(['simulation_', paramSuffix]) = X_m;
    end
end

%% Simulation in the future of the pandemic spread using the model - Prophet estimated parameters

X_prh = zeros(m_states, stepsaheadVec(1));

for jj = 1:length(FitHorizon)
    
    % % starting condition Data
    % X0 = [  y_meas(1,FitHorizon(jj));... S starting cond.
    %         y_meas(2,FitHorizon(jj));... I starting cond.
    %         y_meas(3,FitHorizon(jj));... D starting cond.
    %         y_meas(4,FitHorizon(jj));... T starting cond.
    %         y_meas(5,FitHorizon(jj));... H starting cond.
    %         y_meas(6,FitHorizon(jj))]; % E starting cond.

    % starting condition MHE
    X0 = [  results.sts.S(FitHorizon(jj));... S starting cond.
            results.sts.I(FitHorizon(jj));... I starting cond.
            results.sts.D(FitHorizon(jj));... D starting cond.
            results.sts.T(FitHorizon(jj));... T starting cond.
            results.sts.H(FitHorizon(jj));... H starting cond.
            results.sts.E(FitHorizon(jj))]; % E starting cond.
    
    fieldName = sprintf('wk%d', jj-1);
    ProphParams = getfield(ProphFore, fieldName);

    for ii = 1:length(paramSets)
        X_prh(:,1) = X0;

        for k=1:stepsaheadVec(1)-1
            paramSuffix = paramSets{ii};

            alpha = ProphParams.(['alpha_Prophet_', paramSuffix])(k);
            gamma = ProphParams.(['gamma_Prophet_', paramSuffix])(k);
            delta = ProphParams.(['delta_Prophet_', paramSuffix])(k);
            sigma = ProphParams.(['sigma_Prophet_', paramSuffix])(k);
            tau   = ProphParams.(['tau_Prophet_', paramSuffix])(k);

            eqns2 = [   -m(1) * (alpha * m(2));...
                     m(1) * (alpha * m(2)) - (gamma+lambda) * m(2);...
                     m(2) * gamma - m(3) * (lambda + delta);...
                     delta * m(3) - (tau + sigma)*m(4);...
                     lambda * m(3) + m(4) * sigma + lambda * m(2);...
                     tau * m(4) ];

        f2 = casadi.Function('f2', {m}, {eqns2});

        % Runge-Kutta 4 integration
        k1 = f2(X_prh(:,k));
        k2 = f2(X_prh(:,k) + Ts/2*k1);
        k3 = f2(X_prh(:,k) + Ts/2*k2);
        k4 = f2(X_prh(:,k) + Ts*k3);
        Xprh_plus = Ts/6 * (k1 + 2*k2 + 2*k3 + k4);
        X_prh(:, k+1) = X_prh(:, k) + full(Xprh_plus);
    end
 
        ProphetSim.(fieldName).(['simulation_', paramSuffix]) = X_prh;
    end
end

%% Plotting of previous MHE + Current forecast + real data (Constant Param)
% Infected population 
color1 = [0, 0.4470, 0.7410]; % Blue color
color2 = [0.8500, 0.3250, 0.0980]; % Red color

figure(21)
for jj = 1:length(FitHorizon)
    plt = plot(date( (N_mhe):(N_mhe+FitHorizon(jj)-1) ), results.sts.I(1:FitHorizon(jj)), 's-', 'LineWidth', 1, 'MarkerSize', 5, 'Color', color1);
    plt.MarkerFaceColor = plt.Color;
    hold on
    plot(date(1:FitHorizon(jj)-1), I_data(1:FitHorizon(jj)-1), 'LineWidth', 1.5, 'LineStyle', "-.", 'Color', color2);
    fieldName = sprintf('ForeWk%d_steps21', jj-1); 
    plot(date(FitHorizon(jj):FitHorizon(jj)+stepsaheadVec(1)-1), Xsim.(fieldName)(2,:), 'LineWidth', 2); 
    drawnow;
    pause(1); % Pause for 0.5 seconds
    xlim([date(N_mhe), date(FitHorizon(end)+N_mhe)])
end
plot(date(N_mhe:FitHorizon(end)+N_mhe), I_data(N_mhe:FitHorizon(end)+N_mhe), 'LineWidth', 1.5, 'LineStyle', "-.", 'Color', color2);
hold on 
title('\textbf{I} - Contant parameters','Interpreter','latex')
lgd = legend('MHE fit', 'Real Data', 'Constant \textit{p} Forecast','Interpreter','latex');
lgd.FontSize = 18;

% Infected population  - Holt Prediction
figure(22)
for jj = 1:length(FitHorizon)
    plt = plot(date( (N_mhe):(N_mhe+FitHorizon(jj)-1) ), results.sts.I(1:FitHorizon(jj)), 's-', 'LineWidth', 1, 'MarkerSize', 5, 'Color', color1);
    plt.MarkerFaceColor = plt.Color;
    hold on
    plot(date(1:FitHorizon(jj)-1), I_data(1:FitHorizon(jj)-1), 'LineWidth', 1.5, 'LineStyle', "-.", 'Color', color2);
    fieldName = sprintf('wk%d', jj-1); 
    Xholt_plot=plot(date(FitHorizon(jj):FitHorizon(jj)+stepsaheadVec(1)-1), HoltWintSim.(fieldName).simulation_yhat(2,:), 'LineWidth', 2); 
    Xholt_color = Xholt_plot.Color;
    fill([date(FitHorizon(jj):FitHorizon(jj)+stepsaheadVec(1)-1) flip(date(FitHorizon(jj):FitHorizon(jj)+stepsaheadVec(1)-1))],...
         [HoltWintSim.(fieldName).simulation_yhat_upper(2,:) flip(HoltWintSim.(fieldName).simulation_yhat_lower(2,:))], Xholt_color,'FaceAlpha',.3,'EdgeColor', 'none','HandleVisibility', 'off')
    drawnow;
    pause(1); % Pause for 0.5 seconds
    xlim([date(N_mhe), date(FitHorizon(end)+N_mhe)])
end
hold on 
title('\textbf{I} - Holt Prediction','Interpreter','latex')
lgd = legend('MHE fit', 'Real Data', 'Holt \textit{p} Forecast','Interpreter','latex');
lgd.FontSize = 18;

% Infected population  - Prophet Prediction
figure(23)
for jj = 1:length(FitHorizon)
    plt = plot(date(N_mhe:N_mhe+FitHorizon(jj)-1), results.sts.I(1:FitHorizon(jj)), 's-', 'LineWidth', 1, 'MarkerSize', 5, 'Color', color1);
    plt.MarkerFaceColor = plt.Color;
    hold on
    plot(date(N_mhe:FitHorizon(jj)), I_data(N_mhe:FitHorizon(jj)), 'LineWidth', 1.5, 'LineStyle', "-.", 'Color', color2);
    fieldName = sprintf('wk%d', jj-1); 
    X_prh_plot = plot(date(FitHorizon(jj):FitHorizon(jj)+stepsaheadVec(1)-1), ProphetSim.(fieldName).simulation_yhat(2,:), 'LineWidth', 2);
    X_prh_color = X_prh_plot.Color;
    fill([date(FitHorizon(jj):FitHorizon(jj)+stepsaheadVec(1)-1) flip(date(FitHorizon(jj):FitHorizon(jj)+stepsaheadVec(1)-1))],...
         [ProphetSim.(fieldName).simulation_yhat_upper(2,:) flip(ProphetSim.(fieldName).simulation_yhat_lower(2,:))], X_prh_color,'FaceAlpha',.3,'EdgeColor', 'none','HandleVisibility', 'off')
    drawnow;
    pause(1); % Pause for 0.5 seconds
    xlim([date(N_mhe), date(FitHorizon(end)+N_mhe)])
end
title('\textbf{I} - Prophet predictions','Interpreter','latex')
lgd = legend('MHE fit', 'Real Data', 'Holm \textit{p} Forecast week $n$', 'Interpreter','latex');
lgd.FontSize = 18;

%% Detected population 
figure(31)
for jj = 1:length(FitHorizon)
    plt = plot(date(N_mhe:FitHorizon(jj)), results.sts.D(N_mhe:FitHorizon(jj)), 's-', 'LineWidth', 1, 'MarkerSize', 5, 'Color', color1);
    plt.MarkerFaceColor = plt.Color;
    hold on
    plot(date(N_mhe:FitHorizon(jj)), D_data(N_mhe:FitHorizon(jj)), 'LineWidth', 1.5, 'LineStyle', "-.", 'Color', color2);
    hold on
    fieldName = sprintf('ForeWk%d_steps21', jj-1); 
    plot(date(FitHorizon(jj):FitHorizon(jj)+stepsaheadVec(1)-1), Xsim.(fieldName)(3,:), 'LineWidth', 2); 
    drawnow;
    pause(1); % Pause for 0.5 seconds
    xlim([date(N_mhe), date(FitHorizon(end)+N_mhe)])
end
plot(date(N_mhe:FitHorizon(end)+N_mhe), D_data(N_mhe:FitHorizon(end)+N_mhe), 'LineWidth', 1.5, 'LineStyle', "-.", 'Color', color2);
hold on 
title('\textbf{D} - Contant parameters','Interpreter','latex')
lgd = legend('MHE fit', 'Real Data', 'Constant \textit{p} Forecast','Interpreter','latex');
lgd.FontSize = 18;

% Detected population  - Holt Prediction
figure(32)
for jj = 1:length(FitHorizon)
    plt = plot(date(N_mhe:FitHorizon(jj)), results.sts.D(N_mhe:FitHorizon(jj)), 's-', 'LineWidth', 1, 'MarkerSize', 5, 'Color', color1);
    plt.MarkerFaceColor = plt.Color;
    hold on
    plot(date(N_mhe:FitHorizon(jj)), D_data(N_mhe:FitHorizon(jj)), 'LineWidth', 1.5, 'LineStyle', "-.", 'Color', color2);
    hold on
    fieldName = sprintf('wk%d', jj-1); 
    Xholt_plot=plot(date(FitHorizon(jj):FitHorizon(jj)+stepsaheadVec(1)-1), HoltWintSim.(fieldName).simulation_yhat(3,:), 'LineWidth', 2); 
    hold on
    Xholt_color = Xholt_plot.Color;
    fill([date(FitHorizon(jj):FitHorizon(jj)+stepsaheadVec(1)-1) flip(date(FitHorizon(jj):FitHorizon(jj)+stepsaheadVec(1)-1))],...
         [HoltWintSim.(fieldName).simulation_yhat_upper(3,:) flip(HoltWintSim.(fieldName).simulation_yhat_lower(3,:))], Xholt_color,'FaceAlpha',.3,'EdgeColor', 'none','HandleVisibility', 'off')
    drawnow;
    pause(1); % Pause for 0.5 seconds
    xlim([date(N_mhe), date(FitHorizon(end)+N_mhe)])
end
plot(date(N_mhe:FitHorizon(end)+N_mhe), D_data(N_mhe:FitHorizon(end)+N_mhe), 'LineWidth', 1.5, 'LineStyle', "-.", 'Color', color2);
hold on 
title('\textbf{D} - Holt Prediction','Interpreter','latex')
lgd = legend('MHE fit', 'Real Data', 'Holt \textit{p} Forecast','Interpreter','latex');
lgd.FontSize = 18;

% Detected population - Winters Prediction
figure(33)
for jj = 1:length(FitHorizon)
    plt = plot(date(N_mhe:FitHorizon(jj)), results.sts.D(N_mhe:FitHorizon(jj)), 's-', 'LineWidth', 1, 'MarkerSize', 5, 'Color', color1);
    plt.MarkerFaceColor = plt.Color;
    hold on
    plot(date(N_mhe:FitHorizon(jj)), D_data(N_mhe:FitHorizon(jj)), 'LineWidth', 1.5, 'LineStyle', "-.", 'Color', color2);
    hold on 
    fieldName = sprintf('wk%d', jj-1); 
    X_prh_plot = plot(date(FitHorizon(jj):FitHorizon(jj)+stepsaheadVec(1)-1), ProphetSim.(fieldName).simulation_yhat(3,:), 'LineWidth', 2);
    hold on
    X_prh_color = X_prh_plot.Color;
    fill([date(FitHorizon(jj):FitHorizon(jj)+stepsaheadVec(1)-1) flip(date(FitHorizon(jj):FitHorizon(jj)+stepsaheadVec(1)-1))],...
         [ProphetSim.(fieldName).simulation_yhat_upper(3,:) flip(ProphetSim.(fieldName).simulation_yhat_lower(3,:))], X_prh_color,'FaceAlpha',.3,'EdgeColor', 'none','HandleVisibility', 'off')
    drawnow;
    pause(1); % Pause for 0.5 seconds
    xlim([date(N_mhe), date(FitHorizon(end)+N_mhe)])
end
plot(date(N_mhe:FitHorizon(end)+N_mhe), D_data(N_mhe:FitHorizon(end)+N_mhe), 'LineWidth', 1.5, 'LineStyle', "-.", 'Color', color2);
hold on 
title('\textbf{D} - Prophet predictions','Interpreter','latex')
lgd = legend('MHE fit', 'Real Data', 'Holm \textit{p} Forecast week $n$', 'Interpreter','latex');
lgd.FontSize = 18;

%% Hospitalised population 
figure(41)
for jj = 1:length(FitHorizon)
    plt = plot(date(N_mhe:FitHorizon(jj)), results.sts.T(N_mhe:FitHorizon(jj)), 's-', 'LineWidth', 1, 'MarkerSize', 5, 'Color', color1);
    plt.MarkerFaceColor = plt.Color;
    hold on
    plot(date(N_mhe:FitHorizon(jj)), T_data(N_mhe:FitHorizon(jj)), 'LineWidth', 1.5, 'LineStyle', "-.", 'Color', color2);
    fieldName = sprintf('ForeWk%d_steps21', jj-1); 
    plot(date(FitHorizon(jj):FitHorizon(jj)+stepsaheadVec(1)-1), Xsim.(fieldName)(4,:), 'LineWidth', 2); 
    drawnow;
    pause(1); % Pause for 0.5 seconds
    xlim([date(N_mhe), date(FitHorizon(end)+N_mhe)])
end
title('\textbf{T} - Contant parameters','Interpreter','latex')
lgd = legend('MHE fit', 'Real Data', 'Constant \textit{p} Forecast','Interpreter','latex');
lgd.FontSize = 18;

% Hosp population  - Holt Prediction
figure(42)
for jj = 1:length(FitHorizon)
    plt = plot(date(N_mhe:FitHorizon(jj)), results.sts.T(N_mhe:FitHorizon(jj)), 's-', 'LineWidth', 1, 'MarkerSize', 5, 'Color', color1);
    plt.MarkerFaceColor = plt.Color;
    hold on
    plot(date(N_mhe:FitHorizon(jj)), T_data(N_mhe:FitHorizon(jj)), 'LineWidth', 1.5, 'LineStyle', "-.", 'Color', color2);
    fieldName = sprintf('wk%d', jj-1); 
    Xholt_plot=plot(date(FitHorizon(jj):FitHorizon(jj)+stepsaheadVec(1)-1), HoltWintSim.(fieldName).simulation_yhat(4,:), 'LineWidth', 2); 
    Xholt_color = Xholt_plot.Color;
    fill([date(FitHorizon(jj):FitHorizon(jj)+stepsaheadVec(1)-1) flip(date(FitHorizon(jj):FitHorizon(jj)+stepsaheadVec(1)-1))],...
         [HoltWintSim.(fieldName).simulation_yhat_upper(4,:) flip(HoltWintSim.(fieldName).simulation_yhat_lower(4,:))], Xholt_color,'FaceAlpha',.3,'EdgeColor', 'none','HandleVisibility', 'off')
    drawnow;
    pause(1); % Pause for 0.5 seconds
    xlim([date(N_mhe), date(FitHorizon(end)+N_mhe)])
end
title('\textbf{T} - Holt Prediction','Interpreter','latex')
lgd = legend('MHE fit', 'Real Data', 'Holt \textit{p} Forecast','Interpreter','latex');
lgd.FontSize = 18;


% Hosp population - Winters Prediction
figure(43)
for jj = 1:length(FitHorizon)
    plt = plot(date(N_mhe:FitHorizon(jj)), results.sts.T(N_mhe:FitHorizon(jj)), 's-', 'LineWidth', 1, 'MarkerSize', 5, 'Color', color1);
    plt.MarkerFaceColor = plt.Color;
    hold on
    plot(date(N_mhe:FitHorizon(jj)), T_data(N_mhe:FitHorizon(jj)), 'LineWidth', 1.5, 'LineStyle', "-.", 'Color', color2);
    fieldName = sprintf('wk%d', jj-1); 
    X_prh_plot = plot(date(FitHorizon(jj):FitHorizon(jj)+stepsaheadVec(1)-1), ProphetSim.(fieldName).simulation_yhat(4,:), 'LineWidth', 2);
    X_prh_color = X_prh_plot.Color;
    fill([date(FitHorizon(jj):FitHorizon(jj)+stepsaheadVec(1)-1) flip(date(FitHorizon(jj):FitHorizon(jj)+stepsaheadVec(1)-1))],...
         [ProphetSim.(fieldName).simulation_yhat_upper(4,:) flip(ProphetSim.(fieldName).simulation_yhat_lower(4,:))], X_prh_color,'FaceAlpha',.3,'EdgeColor', 'none','HandleVisibility', 'off')
    drawnow;
    pause(1); % Pause for 0.5 seconds
    xlim([date(N_mhe), date(FitHorizon(end)+N_mhe)])
end
title('\textbf{T} - Prophet predictions','Interpreter','latex')
lgd = legend('MHE fit', 'Real Data', 'Holm \textit{p} Forecast week $n$', 'Interpreter','latex');
lgd.FontSize = 18;

%%
close all
figure()
plot(date(1:21),HoltWintSim.wk0.simulation_yhat(4,:))
hold on 
plot(date(1:21),ProphetSim.wk0.simulation_yhat(4,:))
hold on
plot(date(1:21),T_data(168:168+20))