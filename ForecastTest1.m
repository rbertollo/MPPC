%% Test of forecasting with data as initial conditions
clc
clear all 
close all

load("/Users/marcodelloro/Desktop/MPPC/results.mat")

load('/Users/marcodelloro/Desktop/Thesis/MSc-Thesis-TUe/Data_Collection/Italian Dataset/SIDTTHE_data_DEF.mat');
set(0,'DefaultFigureWindowStyle','docked');

dataset = readtable('/Users/marcodelloro/Desktop/Thesis/MSc-Thesis-TUe/Data_Collection/Italian Dataset/Italy_complete_dataset.xlsx');
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


%% Sim with x0 as Real Data

% for i = 1:size(results.par,1)
% 
%     X0 = [  S_data(i);... S starting cond.
%             I_data(i);... I starting cond.
%             D_data(i);... D starting cond.
%             T_data(i);... T starting cond.
%             H_data(i);... H starting cond.
%             E_data(i)           ];
% 
%      alpha = results.par.alpha(i);
%      gamma = results.par.gamma(i);
%      delta = results.par.delta(i);
%      sigma = results.par.sigma(i);
%      tau   = results.par.tau(i);
% 
%      eqns2 = [ -x(1) * (alpha * x(2));...
%                 x(1) * (alpha * x(2)) - (gamma+lambda) * x(2);...
%                 x(2) * gamma - x(3) * (lambda + delta);...
%                 delta * x(3) - ( tau + sigma )*x(4);...
%                 lambda * x(3) + x(4) * sigma + lambda * x(2);...
%                 tau * x(4)     ];
% 
%     f = casadi.Function('f2', {x}, {eqns2});
% 
%     X_sim(:,1) = X0;
% 
%     for k=1:N_mhe-1
% 
%         % Runge-Kutta 4 integration
%         k1 = f(X_sim(:,k));
%         k2 = f(X_sim(:,k)+Ts/2*k1);
%         k3 = f(X_sim(:,k)+Ts/2*k2);
%         k4 = f(X_sim(:,k)+Ts*k3);
%         x_plus = Ts/6*(k1+2*k2+2*k3+k4);
%         X_sim(:,k+1)=X_sim(:,k) + full(x_plus); % close the gaps - dynamics constraint
%     end
% 
%     Xsim.dyn2 = [Xsim.dyn2, X_sim(:,end)];
% end

%% Plotting of results

% % S 
% figure(1)
% plot(date(1:358),results.stszero.S, 'o') % Initial state
% hold on 
% plot(date(N_mhe:358+N_mhe-1),Xsim.dyn(1,:), LineWidth=1.5) % Initial state
% hold on 
% plot(date(N_mhe:358+N_mhe-1),Xsim.dyn2(1,:), LineWidth=1.5) % Initial state
% hold on 
% plot(date(N_mhe:358+N_mhe-1),results.sts.S, LineWidth=1.5)
% hold on 
% plot(date(1:end),S_data, '--')
% lgd = legend('$\hat{X_0}$', 'Simulated with MHE', 'Simulated with Data', 'MHE results', 'Real Data', 'Interpreter', 'Latex');
% lgd.FontSize = 18;
% xlim([date(1) date(end)])
% 
% % I
% figure(2)
% plot(date(1:358),results.stszero.I, 'o')
% hold on 
% plot(date(N_mhe:358+N_mhe-1),Xsim.dyn(2,:), LineWidth=1.5) % Initial state
% hold on 
% plot(date(N_mhe:358+N_mhe-1),Xsim.dyn2(2,:), LineWidth=1.5) % Initial state
% hold on 
% plot(date(N_mhe:358+N_mhe-1),results.sts.I, LineWidth=1.5)
% hold on 
% plot(date(1:end),I_data, '--')
% lgd = legend('$\hat{X_0}$', 'Simulated with MHE', 'Simulated with Data', 'MHE results', 'Real Data', 'Interpreter', 'Latex');
% lgd.FontSize = 18;
% xlim([date(1) date(end)])
% 
% 
% % D
% figure(3)
% plot(date(1:358),results.stszero.D, 'o')
% hold on 
% plot(date(N_mhe:358+N_mhe-1),Xsim.dyn(3,:), LineWidth=1.5) % Initial state
% hold on 
% plot(date(N_mhe:358+N_mhe-1),Xsim.dyn2(3,:), LineWidth=1.5) % Initial state
% hold on 
% plot(date(N_mhe:358+N_mhe-1),results.sts.D, LineWidth=1.5)
% hold on 
% plot(date(1:end),D_data, '--')
% lgd = legend('$\hat{X_0}$', 'Simulated with MHE', 'Simulated with Data', 'MHE results', 'Real Data', 'Interpreter', 'Latex');
% lgd.FontSize = 18;
% xlim([date(1) date(end)])
% 
% 
% % T
% figure(4)
% plot(date(1:358),results.stszero.T, 'o')
% hold on 
% plot(date(N_mhe:358+N_mhe-1),Xsim.dyn(4,:), LineWidth=1.5) % Initial state
% hold on 
% plot(date(N_mhe:358+N_mhe-1),Xsim.dyn2(4,:), LineWidth=1.5) % Initial state
% hold on 
% plot(date(N_mhe:358+N_mhe-1),results.sts.T, LineWidth=1.5)
% hold on 
% plot(date(1:end),T_data, '--')
% lgd = legend('$\hat{X_0}$', 'Simulated with MHE', 'Simulated with Data', 'MHE results', 'Real Data', 'Interpreter', 'Latex');
% lgd.FontSize = 18;
% xlim([date(1) date(end)])
% 
% 
% % H 
% figure(5)
% plot(date(1:358),results.stszero.H, 'o')
% hold on 
% plot(date(N_mhe:358+N_mhe-1),Xsim.dyn(5,:), LineWidth=1.5) % Initial state
% hold on
% plot(date(N_mhe:358+N_mhe-1),Xsim.dyn2(5,:), LineWidth=1.5) % Initial state
% hold on 
% plot(date(N_mhe:358+N_mhe-1),results.sts.H, LineWidth=1.5)
% hold on 
% plot(date(1:end),H_dataAug, '--')
% lgd = legend('$\hat{X_0}$', 'Simulated with MHE', 'Simulated with Data', 'MHE results', 'Real Data', 'Interpreter', 'Latex');
% lgd.FontSize = 18;
% xlim([date(1) date(end)])


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


for ii = 1:length(paramSets)
    
    paramSuffix = paramSets{ii};
    Xsim.Fore.(sprintf('sim_%s', paramSuffix)) = []

    for j = 1:N_mhe
    
      % X0 = [  S_data(StartHorizon+j+weeksIn(1));... S starting cond.
      %         I_data(StartHorizon+j+weeksIn(1));... I starting cond.
      %         D_data(StartHorizon+j+weeksIn(1));... D starting cond.
      %         T_data(StartHorizon+j+weeksIn(1));... T starting cond.
      %         H_data(StartHorizon+j+weeksIn(1));... H starting cond.
      %         E_data(StartHorizon+j+weeksIn(1))           ];

       X0 = [ results.FullStates{1, 84}.S(j);... S starting cond.
              results.FullStates{1, 84}.I(j);... I starting cond.
              results.FullStates{1, 84}.D(j);...... D starting cond.
              results.FullStates{1, 84}.T(j);...... T starting cond.
              % results.FullStates{1, 84}.H(j);...;... H starting cond.
              results.sts.H(StartHorizon+j+weeksIn(1));... H starting cond.
              results.FullStates{1, 84}.E(j)         ];

        
        alpha = results.par.alpha(StartHorizon+j+weeksIn(1));
        gamma = results.par.gamma(StartHorizon+j+weeksIn(1));    
        % alpha = HoltFore.wk0.(['alpha_', paramSuffix])(j);
        % gamma = HoltFore.wk0.(['gamma_', paramSuffix])(j);
        delta = HoltFore.wk0.(['delta_', paramSuffix])(j);
        sigma = HoltFore.wk0.(['sigma_', paramSuffix])(j);
        tau   = HoltFore.wk0.(['tau_', paramSuffix])(j);
     
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
    
        Xsim.Fore.(sprintf('sim_%s', paramSuffix)) = [Xsim.Fore.(sprintf('sim_%s', paramSuffix)), X_sim(:,end)];
        % ODEres1 = X_sim;
    
    end
end

% Xsim.Fore2 = [];
% Xsim.ForeStart2 = [];

% for j = 1:N_mhe
% 
%   X0 = [  S_data(StartHorizon+j+weeksIn(1));... S starting cond.
%           I_data(StartHorizon+j+weeksIn(1));... I starting cond.
%           D_data(StartHorizon+j+weeksIn(1));... D starting cond.
%           T_data(StartHorizon+j+weeksIn(1));... T starting cond.
%           H_data(StartHorizon+j+weeksIn(1));... H starting cond.
%           E_data(StartHorizon+j+weeksIn(1))           ];
% 
% 
%      alpha = results.par.alpha(StartHorizon+j+weeksIn(1));
%      gamma = results.par.gamma(StartHorizon+j+weeksIn(1));
%      delta = results.par.delta(StartHorizon+j+weeksIn(1));
%      sigma = results.par.sigma(StartHorizon+j+weeksIn(1));
%      tau   = results.par.tau(StartHorizon+j+weeksIn(1));
% 
%      eqns2 = [ -x(1) * (alpha * x(2));...
%                 x(1) * (alpha * x(2)) - (gamma+lambda) * x(2);...
%                 x(2) * gamma - x(3) * (lambda + delta);...
%                 delta * x(3) - ( tau + sigma )*x(4);...
%                 lambda * x(3) + x(4) * sigma + lambda * x(2);...
%                 tau * x(4)     ];
% 
%     f = casadi.Function('f2', {x}, {eqns2});
% 
%     X_sim(:,1) = X0;
% 
%     for k=1:N_mhe-1
% 
%         % Runge-Kutta 4 integration
%         k1 = f(X_sim(:,k));
%         k2 = f(X_sim(:,k)+Ts/2*k1);
%         k3 = f(X_sim(:,k)+Ts/2*k2);
%         k4 = f(X_sim(:,k)+Ts*k3);
%         x_plus = Ts/6*(k1+2*k2+2*k3+k4);
%         X_sim(:,k+1)=X_sim(:,k) + full(x_plus); % close the gaps - dynamics constraint
%     end
% 
%     Xsim.Fore2 = [Xsim.Fore2, X_sim(:,end)];
%     ODEres2 = X_sim;
% 
% end


%% Plotting of results

% I - first week
figure(11)
plot(date(N_mhe:358+N_mhe-1),results.sts.I, LineWidth=1.5)
hold on 
plot(date(1:end),I_data, '--')
hold on 
Xholt_plot=plot(date(StartHorizon+N_mhe:StartHorizon++N_mhe+N_mhe-1), Xsim.Fore.sim_yhat(2,:), LineWidth=1.5);
Xholt_color = Xholt_plot.Color;
hold on
fill([date(StartHorizon+N_mhe:StartHorizon++N_mhe+N_mhe-1) flip(date(StartHorizon+N_mhe:StartHorizon++N_mhe+N_mhe-1))],...
     [Xsim.Fore.sim_yhat_upper(2,:) flip(Xsim.Fore.sim_yhat_lower(2,:))], Xholt_color,'FaceAlpha',.3,'EdgeColor', 'none','HandleVisibility', 'off')
hold on
xline(date(StartHorizon+N_mhe),":",'HandleVisibility', 'off')
lgd = legend('MHE results', 'Real Data', '21 days Future Forecast','Simulation From Data','21 days Future Forecast sim param','Interpreter', 'Latex');
lgd.FontSize = 18;
xlim([date(1) date(end)])


% D - first week
figure(12)
plot(date(N_mhe:358+N_mhe-1),results.sts.D, LineWidth=1.5)
hold on 
plot(date(1:end),D_data, '--')
hold on 
Xholt_plot=plot(date(StartHorizon+N_mhe:StartHorizon++N_mhe+N_mhe-1), Xsim.Fore.sim_yhat(3,:), LineWidth=1.5);
Xholt_color = Xholt_plot.Color;
hold on
fill([date(StartHorizon+N_mhe:StartHorizon++N_mhe+N_mhe-1) flip(date(StartHorizon+N_mhe:StartHorizon++N_mhe+N_mhe-1))],...
     [Xsim.Fore.sim_yhat_upper(3,:) flip(Xsim.Fore.sim_yhat_lower(3,:))], Xholt_color,'FaceAlpha',.3,'EdgeColor', 'none','HandleVisibility', 'off')
hold on
plot(date(N_mhe:358+N_mhe-1),Xsim.dyn2(3,:), LineWidth=1.5) % New fitted model
hold on 
plot(date(StartHorizon+N_mhe:StartHorizon+N_mhe+N_mhe-1), Xsim.Fore2(3,:), LineWidth=2 )
hold on 
xline(date(StartHorizon+N_mhe),":",'HandleVisibility', 'off')
lgd = legend('MHE results', 'Real Data', '21 days Future Forecast','Simulation From Data','21 days Future Forecast sim param','Interpreter', 'Latex');
lgd.FontSize = 18;
xlim([date(1) date(end)])


% T - first week
figure(13)
plot(date(N_mhe:358+N_mhe-1),results.sts.T, LineWidth=1.5)
hold on 
plot(date(1:end),T_data, '--')
hold on 
Xholt_plot=plot(date(StartHorizon+N_mhe:StartHorizon++N_mhe+N_mhe-1), Xsim.Fore.sim_yhat(4,:), LineWidth=1.5);
Xholt_color = Xholt_plot.Color;
hold on
fill([date(StartHorizon+N_mhe:StartHorizon+N_mhe+N_mhe-1) flip(date(StartHorizon+N_mhe:StartHorizon+N_mhe+N_mhe-1))],...
     [Xsim.Fore.sim_yhat_upper(4,:) flip(Xsim.Fore.sim_yhat_lower(4,:))], Xholt_color,'FaceAlpha',.3,'EdgeColor', 'none','HandleVisibility', 'off')
hold on
plot(date(N_mhe:358+N_mhe-1),Xsim.dyn2(4,:), LineWidth=1.5) % New fitted model
hold on 
plot(date(StartHorizon+N_mhe:StartHorizon+N_mhe+N_mhe-1), Xsim.Fore2(4,:), LineWidth=2 )
hold on 
xline(date(StartHorizon+N_mhe),":",'HandleVisibility', 'off')
lgd = legend('MHE results', 'Real Data', '21 days Future Forecast','Simulation From Data','21 days Future Forecast sim param','Interpreter', 'Latex');
lgd.FontSize = 18;
xlim([date(1) date(end)])


% H - first week
figure(14)
plot(date(N_mhe:358+N_mhe-1),results.sts.H, LineWidth=1.5)
hold on 
plot(date(1:end),H_data, '--')
hold on 
Xholt_plot=plot(date(StartHorizon+N_mhe:StartHorizon++N_mhe+N_mhe-1), Xsim.Fore.sim_yhat(5,:), LineWidth=1.5);
Xholt_color = Xholt_plot.Color;
hold on
fill([date(StartHorizon+N_mhe:StartHorizon+N_mhe+N_mhe-1) flip(date(StartHorizon+N_mhe:StartHorizon+N_mhe+N_mhe-1))],...
     [Xsim.Fore.sim_yhat_upper(5,:) flip(Xsim.Fore.sim_yhat_lower(5,:))], Xholt_color,'FaceAlpha',.3,'EdgeColor', 'none','HandleVisibility', 'off')
hold on
xline(date(StartHorizon+N_mhe),":",'HandleVisibility', 'off')
lgd = legend('MHE results', 'Real Data', '21 days Future Forecast','Simulation From Data','21 days Future Forecast sim param','Interpreter', 'Latex');
lgd.FontSize = 18;
xlim([date(1) date(end)])


% E - first week
figure(15)
plot(date(N_mhe:358+N_mhe-1),results.sts.E, LineWidth=1.5)
hold on 
plot(date(1:end),E_data, '--')
hold on 
Xholt_plot=plot(date(StartHorizon+N_mhe:StartHorizon++N_mhe+N_mhe-1), Xsim.Fore.sim_yhat(6,:), LineWidth=1.5);
Xholt_color = Xholt_plot.Color;
hold on
fill([date(StartHorizon+N_mhe:StartHorizon+N_mhe+N_mhe-1) flip(date(StartHorizon+N_mhe:StartHorizon+N_mhe+N_mhe-1))],...
     [Xsim.Fore.sim_yhat_upper(6,:) flip(Xsim.Fore.sim_yhat_lower(6,:))], Xholt_color,'FaceAlpha',.3,'EdgeColor', 'none','HandleVisibility', 'off')
hold on
plot(date(N_mhe:358+N_mhe-1),Xsim.dyn2(6,:), LineWidth=1.5) % New fitted model
hold on 
plot(date(StartHorizon+N_mhe:StartHorizon+N_mhe+N_mhe-1), Xsim.Fore2(6,:), LineWidth=2 )
hold on 
xline(date(StartHorizon+N_mhe),":",'HandleVisibility', 'off')
lgd = legend('MHE results', 'Real Data', '21 days Future Forecast','Simulation From Data','21 days Future Forecast sim param','Interpreter', 'Latex');
lgd.FontSize = 18;
xlim([date(1) date(end)])


%% Initial values check / Trial with different starting conditions

dynTry.Fore = [];
dynTry.MHE = [];

for j = 0:N_mhe-1

X0 = [ results.stszero.S(StartHorizon+j+weeksIn(1));... S starting cond.
            results.stszero.I(StartHorizon+j+weeksIn(1));... I starting cond.
            results.stszero.D(StartHorizon+j+weeksIn(1));... D starting cond.
            results.stszero.T(StartHorizon+j+weeksIn(1));... T starting cond.
            results.stszero.H(StartHorizon+j+weeksIn(1));... H starting cond.
            results.stszero.E(StartHorizon+j+weeksIn(1))          ];

            alpha = results.par.alpha(StartHorizon+j+weeksIn(1));
            gamma = results.par.gamma(StartHorizon+j+weeksIn(1));
            delta = results.par.delta(StartHorizon+j+weeksIn(1));
            sigma = results.par.sigma(StartHorizon+j+weeksIn(1));
            tau   = results.par.tau(StartHorizon+j+weeksIn(1));

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

    dynTry.MHE = [dynTry.MHE, X_sim(:,end)];
    ODEres2 = X_sim;
 
end

% Dynamics on the last values from MHE w/ "Real" coefficients

estStatesFore = results.FullStates{1, 84};

for j = 0:N_mhe-1

X0 = [ estStatesFore.S(j+1);... S starting cond.
       estStatesFore.I(j+1);... I starting cond.
       estStatesFore.D(j+1);... D starting cond.
       estStatesFore.T(j+1);... T starting cond.
       estStatesFore.H((j+1));... H starting cond.
       estStatesFore.E((j+1))          ];

    alpha = results.par.alpha(StartHorizon+j+weeksIn(1));
    gamma = results.par.gamma(StartHorizon+j+weeksIn(1));
    delta = results.par.delta(StartHorizon+j+weeksIn(1));
    sigma = results.par.sigma(StartHorizon+j+weeksIn(1));
    tau   = results.par.tau(StartHorizon+j+weeksIn(1));

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

% Dynamics on the last values from MHE w/ Forecasted coefficients
dynTry.Fore2 = [];
for j = 0:N_mhe-1

X0 = [ estStatesFore.S(j+1);... S starting cond.
       estStatesFore.I(j+1);... I starting cond.
       estStatesFore.D(j+1);... D starting cond.
       estStatesFore.T(j+1);... T starting cond.
       estStatesFore.H((j+1));... H starting cond.
       estStatesFore.E((j+1))          ];

    alpha = HoltFore.wk0.alpha_yhat(j+1);
    gamma = HoltFore.wk0.gamma_yhat(j+1);
    delta = HoltFore.wk0.delta_yhat(j+1);
    sigma = HoltFore.wk0.sigma_yhat(j+1);
    tau   = HoltFore.wk0.tau_yhat(j+1);

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

    dynTry.Fore2 = [dynTry.Fore2, X_sim(:,end)];
    ODEres2 = X_sim;
 
end

figure(50)
plot( date(StartHorizon:StartHorizon+N_mhe-1), dynTry.MHE(1,:), LineWidth=2 );
hold on
plot( date(StartHorizon:StartHorizon+N_mhe-1), dynTry.Fore(1,:), LineWidth=2 );
hold on
plot( date(StartHorizon:StartHorizon+N_mhe-1), dynTry.Fore2(1,:), LineWidth=2 );
hold on
plot( date(StartHorizon:StartHorizon+N_mhe-1), results.sts.S(StartHorizon:StartHorizon+N_mhe-1),LineWidth=2);
lgd = legend('Simulated from $X_0$ MHE', 'Simulated With Full Results REAL $p$', 'Simulated With Full Results FORECAST $p$', 'MHE Results','Interpreter', 'Latex');
lgd.FontSize = 18;


figure(51)
plot( date(StartHorizon:StartHorizon+N_mhe-1), dynTry.MHE(2,:), LineWidth=2 );
hold on
plot( date(StartHorizon:StartHorizon+N_mhe-1), dynTry.Fore(2,:), LineWidth=2 );
hold on
plot( date(StartHorizon:StartHorizon+N_mhe-1), dynTry.Fore2(2,:), LineWidth=2 );
hold on
plot( date(StartHorizon:StartHorizon+N_mhe-1), results.sts.I(StartHorizon:StartHorizon+N_mhe-1),LineWidth=2);
lgd = legend('Simulated from $X_0$ MHE', 'Simulated With Full Results REAL $p$', 'Simulated With Full Results FORECAST $p$', 'MHE Results','Interpreter', 'Latex');
lgd.FontSize = 18;


figure(52)
plot( date(StartHorizon:StartHorizon+N_mhe-1), dynTry.MHE(3,:), LineWidth=2 );
hold on
plot( date(StartHorizon:StartHorizon+N_mhe-1), dynTry.Fore(3,:), LineWidth=2 );
hold on
plot( date(StartHorizon:StartHorizon+N_mhe-1), dynTry.Fore2(3,:), LineWidth=2 );
hold on
plot( date(StartHorizon:StartHorizon+N_mhe-1), results.sts.D(StartHorizon:StartHorizon+N_mhe-1),LineWidth=2);
lgd = legend('Simulated from $X_0$ MHE', 'Simulated With Full Results REAL $p$', 'Simulated With Full Results FORECAST $p$', 'MHE Results','Interpreter', 'Latex');
lgd.FontSize = 18;

figure(53)
plot( date(StartHorizon:StartHorizon+N_mhe-1), dynTry.MHE(4,:), LineWidth=2 );
hold on
plot( date(StartHorizon:StartHorizon+N_mhe-1), dynTry.Fore(4,:), LineWidth=2 );
hold on
plot( date(StartHorizon:StartHorizon+N_mhe-1), dynTry.Fore2(4,:), LineWidth=2 );
hold on
plot( date(StartHorizon:StartHorizon+N_mhe-1), results.sts.T(StartHorizon:StartHorizon+N_mhe-1),LineWidth=2);
lgd = legend('Simulated from $X_0$ MHE', 'Simulated With Full Results REAL $p$', 'Simulated With Full Results FORECAST $p$', 'MHE Results','Interpreter', 'Latex');
lgd.FontSize = 18;

