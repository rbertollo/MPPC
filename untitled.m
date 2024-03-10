%% Script with Forecast of 5 weeks - Last version
clc
clear all 
close all

load("/Users/marcodelloro/Desktop/MPPC/results.mat")
load('/Users/marcodelloro/Desktop/Thesis/MSc-Thesis-TUe/Data_Collection/Italian Dataset/SIDTTHE_data_DEF.mat');
set(0,'DefaultFigureWindowStyle','docked');

dataset = readtable('/Users/marcodelloro/Desktop/Thesis/MSc-Thesis-TUe/Data_Collection/Italian Dataset/Italy_complete_dataset.xlsx');
dataset = dataset(190:588,:);
dataset.data = datetime(dataset.data, 'InputFormat', 'dd-MMM-yyyy HH:mm:ss');

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
lambda = 0.1;

x = [s; i; d; t; h; e];

n_states = length(x);

Ts = 1;
Xsim.dyn = [];
Xsim.dyn2 = [];

%% Sim with x0 as results from MHE

for i = 1:size(results.par,1)

    X0 = [  results.stszero.S(i);... S starting cond.
            results.stszero.I(i);... I starting cond.
            results.stszero.D(i);... D starting cond.
            results.stszero.T(i);... T starting cond.
            results.stszero.H(i);... H starting cond.
            results.stszero.E(i) ];

     alpha = results.par.alpha(i);
     gamma = results.par.gamma(i);
     delta = results.par.delta(i);
     sigma = results.par.sigma(i);
     tau   = results.par.tau(i);

     eqns2 = [ -x(1) * (alpha * x(2));...
                x(1) * (alpha * x(2)) - (gamma+lambda) * x(2);...
                x(2) * gamma - x(3) * (lambda + delta);...
                delta * x(3) - ( tau + sigma )*x(4);...
                lambda * x(3) + x(4) * sigma + lambda * x(2);...
                tau * x(4)     ];

    f = casadi.Function('f2', {x}, {eqns2});
   
    X_sim(:,1) = X0;

    for k=1:N_mhe-1
       
        % Runge-Kutta 4 integration
        k1 = f(X_sim(:,k));
        k2 = f(X_sim(:,k)+Ts/2*k1);
        k3 = f(X_sim(:,k)+Ts/2*k2);
        k4 = f(X_sim(:,k)+Ts*k3);
        x_plus = Ts/6*(k1+2*k2+2*k3+k4);
        X_sim(:,k+1)=X_sim(:,k) + full(x_plus); % close the gaps - dynamics constraint
    end
   
    Xsim.dyn = [Xsim.dyn, X_sim(:,end)];
end


%% Try to simulate in the future 

ForeWk0 = readtable('/Users/marcodelloro/Desktop/Thesis/MSc-Thesis-TUe/Forecasting/TrialPython/Forecasts21d_CI2.xlsx', 'Sheet', 'Horizon_0_days');
ForeWk1 = readtable('/Users/marcodelloro/Desktop/Thesis/MSc-Thesis-TUe/Forecasting/TrialPython/Forecasts21d_CI2.xlsx', 'Sheet', 'Horizon_7_days');
ForeWk2 = readtable('/Users/marcodelloro/Desktop/Thesis/MSc-Thesis-TUe/Forecasting/TrialPython/Forecasts21d_CI2.xlsx', 'Sheet', 'Horizon_14_days');
ForeWk3 = readtable('/Users/marcodelloro/Desktop/Thesis/MSc-Thesis-TUe/Forecasting/TrialPython/Forecasts21d_CI2.xlsx', 'Sheet', 'Horizon_21_days');
ForeWk4 = readtable('/Users/marcodelloro/Desktop/Thesis/MSc-Thesis-TUe/Forecasting/TrialPython/Forecasts21d_CI2.xlsx', 'Sheet', 'Horizon_28_days');
ForeWk5 = readtable('/Users/marcodelloro/Desktop/Thesis/MSc-Thesis-TUe/Forecasting/TrialPython/Forecasts21d_CI2.xlsx', 'Sheet', 'Horizon_35_days');

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

StartHorizon = 4*N_mhe; % we are here
weeksIn = [0, 7, 14, 21, 28, 35];

% Forecast for the next 21 day from where we are w/ Data ()
Xsim.Fore = [];
Xsim.ForeStart = [];

paramSets = {'yhat', 'yhat_lower', 'yhat_upper'};
% Dynamics on the last values from MHE w/ "Real" coefficients

%% Evaluataion of H to see what is ruining the estimation

estStatesFore = results.FullStates{1, 84};
dynTry.Fore = [];

for j = 0:N_mhe-1

X0 = [  estStatesFore.S(j+1);... S starting cond.
        estStatesFore.I(j+1);... I starting cond.
        estStatesFore.D(j+1);... D starting cond.
        estStatesFore.T(j+1);... T starting cond.
        estStatesFore.H((j+1));... H starting cond.
        % estStatesFore.E(j+1)

        % results.stszero.S(StartHorizon+j+weeksIn(1));... S starting cond.
        % results.stszero.I(StartHorizon+j+weeksIn(1));... I starting cond.
        % results.stszero.D(StartHorizon+j+weeksIn(1));... D starting cond.
        % results.stszero.T(StartHorizon+j+weeksIn(1));... D starting cond.
        % estStatesFore.H((j+1));... H starting cond.
        % results.stszero.H(StartHorizon+j+weeksIn(1));... T starting cond.
        results.stszero.E(StartHorizon+weeksIn(1))
        ];

     % alpha = results.par.alpha(StartHorizon+j+weeksIn(1));
     % gamma = results.par.gamma(StartHorizon+j+weeksIn(1));
     % delta = results.par.delta(StartHorizon+j+weeksIn(1));
     % sigma = results.par.sigma(StartHorizon+j+weeksIn(1));
     % tau   = results.par.tau(StartHorizon+j+weeksIn(1));

     alpha = HoltFore.wk0.alpha_yhat(j+1);
     gamma = HoltFore.wk0.gamma_yhat(j+1);
     delta = HoltFore.wk0.delta_yhat(j+1);
     sigma = HoltFore.wk0.sigma_yhat(j+1);
     tau = HoltFore.wk0.tau_yhat(j+1);
   


     eqns2 = [ -x(1) * (alpha * x(2));...
                x(1) * (alpha * x(2)) - (gamma+lambda) * x(2);...
                x(2) * gamma - x(3) * (lambda + delta);...
                delta * x(3) - ( tau + sigma )*x(4);...
                lambda * x(3) + x(4) * sigma + lambda * x(2);...
                tau * x(4)     ];

    f = casadi.Function('f2', {x}, {eqns2});

    X_sim(:,1) = X0;

    for k=1:N_mhe-1

        % Runge-Kutta 4 integration
        k1 = f(X_sim(:,k));
        k2 = f(X_sim(:,k)+Ts/2*k1);
        k3 = f(X_sim(:,k)+Ts/2*k2);
        k4 = f(X_sim(:,k)+Ts*k3);
        x_plus = Ts/6*(k1+2*k2+2*k3+k4);
        X_sim(:,k+1)=X_sim(:,k) + full(x_plus); % close the gaps - dynamics constraint
    end

    dynTry.Fore = [dynTry.Fore, X_sim(:,end)];
    ODEres2 = X_sim;
 
end

% Plot of the H x_0 (starting condition)

 startIdx = N_mhe + StartHorizon + weeksIn(1) - 1;
 endIdx = startIdx + N_mhe - 1;

%% E
figure(14)
plot(date(1:end), E_data, '--') % Plot the real data
hold on
plot(date(N_mhe:358+N_mhe-1), Xsim.dyn(6,:), 'LineWidth', 1.5) % Plot the new fitted model
hold on
plot(date(N_mhe:358+N_mhe-1), results.stszero.E, 'LineWidth', 1.5) % Plot the new fitted model
hold on
Xholt_plot = plot(date(startIdx:endIdx ), dynTry.Fore(6,:), 'LineWidth', 2.5);
Xholt_color = Xholt_plot.Color;
hold on
xline(date(StartHorizon+N_mhe),":",'HandleVisibility', 'off')
lgd = legend('Real Data', 'MHE results', sprintf('%d days Future Forecast', N_mhe), 'Interpreter', 'Latex');
lgd.FontSize = 18;
xlim([date(1) date(end)])
title('E-Expired', 'Interpreter', 'Latex'); 

% Plot of the H x_0 (starting condition)
figure(15)
plot(date(1:N_mhe+1), results.stszero.E(StartHorizon:StartHorizon+N_mhe),'LineWidth', 1.5) 
hold on
plot(date(1:N_mhe),  estStatesFore.E ,'LineWidth', 1.5) 
lgd = legend('Good $x_0$', 'Assumed $x_0$', 'Interpreter','Latex','location','southeast');
lgd.FontSize = 18;


%% Score Function 

% Score on E-deceased population

for pp = 1:1
    startIdx = N_mhe + StartHorizon + weeksIn(pp) - 1;
    startIdxMHE = StartHorizon + weeksIn(pp);

    currentWeekField = sprintf('wk%d', pp-1);
    weekend = [7,14,21];
    
    for ii = 1:3
        % score week by week for holt model
        SE = ( dynTry.Fore.sim_yhat(6,weekend(ii)) - E_data(startIdx + weekend(ii)) ) ./ ( E_data(startIdx + weekend(ii)) -  E_data(startIdx) );
        score.holtData.(currentWeekField).(['week', num2str(weekend(ii))]) = SE_holt;
        score.holtMHE.(currentWeekField).(['week', num2str(weekend(ii))]) = SE_holtMHE;
        score.ProphMHE.(currentWeekField).(['week', num2str(weekend(ii))]) = SE_ProphMHE;

        % baseline score
        SE_baseline = ( 2*E_data(startIdx) - E_data(startIdx - weekend(1)) - E_data(startIdx + weekend(ii)) ) ./ ( E_data(startIdx + weekend(ii)) -  E_data(startIdx) );
        SE_BaseMHE = ( 2*Xsim.dyn(6,startIdxMHE) - Xsim.dyn(6,startIdxMHE - weekend(1)) - Xsim.dyn(6,startIdxMHE + weekend(ii)) ) ./ ( Xsim.dyn(6,startIdxMHE + weekend(ii)) -  Xsim.dyn(6,startIdxMHE) );
        score.baseData.(currentWeekField).(['week', num2str(weekend(ii))]) = SE_baseline;
        score.baseMHE.(currentWeekField).(['week', num2str(weekend(ii))]) = SE_BaseMHE;
    end
end