%% Script with Forecast of 5 weeks - Third time frame
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

%% Try to simulate in the future 
% Load of Forecasts21d_CI4.xlsx (horizon from 7*N_mhe)

ForeWk0 = readtable('/Users/marcodelloro/Desktop/Thesis/MSc-Thesis-TUe/Forecasting/TrialPython/Forecasts21d_CI4.xlsx', 'Sheet', 'Horizon_0_days');
ForeWk1 = readtable('/Users/marcodelloro/Desktop/Thesis/MSc-Thesis-TUe/Forecasting/TrialPython/Forecasts21d_CI4.xlsx', 'Sheet', 'Horizon_7_days');
ForeWk2 = readtable('/Users/marcodelloro/Desktop/Thesis/MSc-Thesis-TUe/Forecasting/TrialPython/Forecasts21d_CI4.xlsx', 'Sheet', 'Horizon_14_days');
ForeWk3 = readtable('/Users/marcodelloro/Desktop/Thesis/MSc-Thesis-TUe/Forecasting/TrialPython/Forecasts21d_CI4.xlsx', 'Sheet', 'Horizon_21_days');
ForeWk4 = readtable('/Users/marcodelloro/Desktop/Thesis/MSc-Thesis-TUe/Forecasting/TrialPython/Forecasts21d_CI4.xlsx', 'Sheet', 'Horizon_28_days');
ForeWk5 = readtable('/Users/marcodelloro/Desktop/Thesis/MSc-Thesis-TUe/Forecasting/TrialPython/Forecasts21d_CI4.xlsx', 'Sheet', 'Horizon_35_days');

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

StartHorizon = 7*N_mhe; % we are here --> NEW HORIZON
weeksIn = [0, 7, 14, 21, 28, 35];

% Forecast for the next 21 day from where we are w/ Data ()
Xsim.Fore = [];
Xsim.ForeStart = [];

paramSets = {'yhat', 'yhat_lower', 'yhat_upper'};
% Dynamics on the last values from MHE w/ "Real" coefficients

%% Constant Parameter forecast
for pp = 1:length(weeksIn)
    
    estStatesFore = results.FullStates{1, StartHorizon+weeksIn(pp)};
    weekField = sprintf('wk%d', pp-1);
        
    Sim.Const.(weekField) = [];
    
        for j = 0:N_mhe-1
    
            X0 = [ estStatesFore.S(j+1);... S starting cond.
                   estStatesFore.I(j+1);... I starting cond.
                   estStatesFore.D(j+1);... D starting cond.
                   estStatesFore.T(j+1);... T starting cond.
                   estStatesFore.H(j+1);... H starting cond.
                   estStatesFore.E(j+1)          ];
                
            alpha = HoltFore.(weekField).alpha_yhat(1);
            gamma = HoltFore.(weekField).gamma_yhat(1);
            delta = HoltFore.(weekField).delta_yhat(1);
            sigma = HoltFore.(weekField).sigma_yhat(1);
            tau   = HoltFore.(weekField).tau_yhat(1);
    
    
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
        
            Sim.Const.(weekField) = [Sim.Const.(weekField), X_sim(:, end)];
        
        end
end

%% Holt Parameter forecast
for pp = 1:length(weeksIn)
    
    estStatesFore = results.FullStates{1, StartHorizon+weeksIn(pp)};
    weekField = sprintf('wk%d', pp-1);
    
    for ii = 1:length(paramSets)
        
        paramSuffix = paramSets{ii};
        Sim.Holt.(weekField).(sprintf('sim_%s', paramSuffix)) = [];
    
        for j = 0:N_mhe-1
    
            X0 = [ estStatesFore.S(j+1);... S starting cond.
                   estStatesFore.I(j+1);... I starting cond.
                   estStatesFore.D(j+1);... D starting cond.
                   estStatesFore.T(j+1);... T starting cond.
                   estStatesFore.H(j+1);... H starting cond.
                   estStatesFore.E(j+1)          ];
                
            alpha = HoltFore.(weekField).(['alpha_', paramSuffix])(j+1);
            gamma = HoltFore.(weekField).(['gamma_', paramSuffix])(j+1);
            delta = HoltFore.(weekField).(['delta_', paramSuffix])(j+1);
            sigma = HoltFore.(weekField).(['sigma_', paramSuffix])(j+1);
            tau   = HoltFore.(weekField).(['tau_', paramSuffix])(j+1);
    
    
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
        
            Sim.Holt.(weekField).(sprintf('sim_%s', paramSuffix)) = [Sim.Holt.(weekField).(sprintf('sim_%s', paramSuffix)), X_sim(:, end)];
        
        end
    end
end

%% Prophet Parameter forecast
for pp = 1:length(weeksIn)
    
    estStatesFore = results.FullStates{1, StartHorizon+weeksIn(pp)};
    weekField = sprintf('wk%d', pp-1);
    
    for ii = 1:length(paramSets)
        
        paramSuffix = paramSets{ii};
        Sim.Proph.(weekField).(sprintf('sim_%s', paramSuffix)) = [];
    
        for j = 0:N_mhe-1
    
            X0 = [ estStatesFore.S(j+1);... S starting cond.
                   estStatesFore.I(j+1);... I starting cond.
                   estStatesFore.D(j+1);... D starting cond.
                   estStatesFore.T(j+1);... T starting cond.
                   estStatesFore.H(j+1);... H starting cond.
                   estStatesFore.E(j+1)          ];
                
            alpha = ProphFore.(weekField).(['alpha_Prophet_', paramSuffix])(j+1);
            gamma = ProphFore.(weekField).(['gamma_Prophet_', paramSuffix])(j+1);
            delta = ProphFore.(weekField).(['delta_Prophet_', paramSuffix])(j+1);
            sigma = ProphFore.(weekField).(['sigma_Prophet_', paramSuffix])(j+1);
            tau   = ProphFore.(weekField).(['tau_Prophet_', paramSuffix])(j+1);
    
    
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
        
            Sim.Proph.(weekField).(sprintf('sim_%s', paramSuffix)) = [Sim.Proph.(weekField).(sprintf('sim_%s', paramSuffix)), X_sim(:, end)];
        
        end
    end
end



%% Plotting of results (Holt)

% I 
figure(11)
plot(date(1:end), I_data, '--') % Plot the real data
hold on
plot(date(N_mhe:358+N_mhe-1), Xsim.dyn(2,:), 'LineWidth', 1.5) % Plot the new fitted model

for pp = 1:length(weeksIn)
    currentWeekField = sprintf('wk%d', pp-1);
    startIdx = N_mhe + StartHorizon + weeksIn(pp) - 1;
    endIdx = startIdx + N_mhe - 1;
    
    Xholt_plot = plot(date(startIdx:endIdx), Sim.Holt.(currentWeekField).sim_yhat(2,:), 'LineWidth', 1.5);
    Xholt_color = Xholt_plot.Color;
    fill([date(startIdx:endIdx) flip(date(startIdx:endIdx))],...
         [Sim.Holt.(currentWeekField).sim_yhat_upper(2,:) flip(Sim.Holt.(currentWeekField).sim_yhat_lower(2,:))],...
         Xholt_color, 'FaceAlpha', .3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    hold on
    Xcost_plt = plot(date(startIdx:endIdx), Sim.Const.(currentWeekField)(2,:), 'LineWidth', 2,'LineStyle','-.');

    xline(date(startIdx), ":", 'HandleVisibility', 'off');
    
    hold on 
end
lgd = legend('Real Data', 'MHE results', sprintf('%d days Future Forecast', N_mhe), 'Interpreter', 'Latex');
lgd.FontSize = 18;
xlim([date(1) date(end)])
title('I-Infected', 'Interpreter', 'Latex');


% D 
figure(12)
plot(date(1:end), D_data, '--') % Plot the real data
hold on
plot(date(N_mhe:358+N_mhe-1), Xsim.dyn(3,:), 'LineWidth', 1.5) % Plot the new fitted model

for pp = 1:length(weeksIn)
    currentWeekField = sprintf('wk%d', pp-1);
    startIdx = N_mhe + StartHorizon + weeksIn(pp) - 1;
    endIdx = startIdx + N_mhe - 1;
    
    Xholt_plot = plot(date(startIdx:endIdx), Sim.Holt.(currentWeekField).sim_yhat(3,:), 'LineWidth', 1.5);
    Xholt_color = Xholt_plot.Color;
    fill([date(startIdx:endIdx) flip(date(startIdx:endIdx))],...
         [Sim.Holt.(currentWeekField).sim_yhat_upper(3,:) flip(Sim.Holt.(currentWeekField).sim_yhat_lower(3,:))],...
         Xholt_color, 'FaceAlpha', .3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    hold on
    Xcost_plt = plot(date(startIdx:endIdx), Sim.Const.(currentWeekField)(3,:), 'LineWidth', 2,'LineStyle','-.');


    xline(date(startIdx), ":", 'HandleVisibility', 'off');
    
    hold on
end
lgd = legend('Real Data', 'MHE results', sprintf('%d days Future Forecast', N_mhe), 'Interpreter', 'Latex');
lgd.FontSize = 18;
xlim([date(1) date(end)])
title('D-Infected', 'Interpreter', 'Latex'); 


% T 
figure(13)
plot(date(1:end), T_data, '--') % Plot the real data
hold on
plot(date(N_mhe:358+N_mhe-1), Xsim.dyn(4,:), 'LineWidth', 1.5) % Plot the new fitted model

for pp = 1:length(weeksIn)
    currentWeekField = sprintf('wk%d', pp-1);
    startIdx = N_mhe + StartHorizon + weeksIn(pp) - 1;
    endIdx = startIdx + N_mhe - 1;
    
    Xholt_plot = plot(date(startIdx:endIdx), Sim.Holt.(currentWeekField).sim_yhat(4,:), 'LineWidth', 1.5);
    Xholt_color = Xholt_plot.Color;
    fill([date(startIdx:endIdx) flip(date(startIdx:endIdx))],...
         [Sim.Holt.(currentWeekField).sim_yhat_upper(4,:) flip(Sim.Holt.(currentWeekField).sim_yhat_lower(4,:))],...
         Xholt_color, 'FaceAlpha', .3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    hold on
    Xcost_plt = plot(date(startIdx:endIdx), Sim.Const.(currentWeekField)(4,:), 'LineWidth', 2,'LineStyle','-.');

    xline(date(startIdx), ":", 'HandleVisibility', 'off');
    
    hold on
end
lgd = legend('Real Data', 'MHE results', sprintf('%d days Future Forecast', N_mhe), 'Interpreter', 'Latex');
lgd.FontSize = 18;
xlim([date(1) date(end)])
title('T-Threatned', 'Interpreter', 'Latex'); 


% H 
figure(14)
plot(date(1:end), H_data, '--') % Plot the real data
hold on
plot(date(N_mhe:358+N_mhe-1), Xsim.dyn(5,:), 'LineWidth', 1.5) % Plot the new fitted model

for pp = 1:length(weeksIn)
    currentWeekField = sprintf('wk%d', pp-1);
    startIdx = N_mhe + StartHorizon + weeksIn(pp) - 1;
    endIdx = startIdx + N_mhe - 1;
    
    Xholt_plot = plot(date(startIdx:endIdx), Sim.Holt.(currentWeekField).sim_yhat(5,:), 'LineWidth', 1.5);
    Xholt_color = Xholt_plot.Color;
    fill([date(startIdx:endIdx) flip(date(startIdx:endIdx))],...
         [Sim.Holt.(currentWeekField).sim_yhat_upper(5,:) flip(Sim.Holt.(currentWeekField).sim_yhat_lower(5,:))],...
         Xholt_color, 'FaceAlpha', .3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    hold on
    Xcost_plt = plot(date(startIdx:endIdx), Sim.Const.(currentWeekField)(5,:), 'LineWidth', 2,'LineStyle','-.');

    xline(date(startIdx), ":", 'HandleVisibility', 'off');
    
    hold on 
end
lgd = legend('Real Data', 'MHE results', sprintf('%d days Future Forecast', N_mhe), 'Interpreter', 'Latex');
lgd.FontSize = 18;
xlim([date(1) date(end)])
title('H-Healed', 'Interpreter', 'Latex'); 

% E - first week
figure(15)
plot(date(1:end), E_data, '--') % Plot the real data
hold on
plot(date(N_mhe:358+N_mhe-1), Xsim.dyn(6,:), 'LineWidth', 1.5) % Plot the new fitted model

for pp = 1:length(weeksIn)
    currentWeekField = sprintf('wk%d', pp-1);
    startIdx = N_mhe + StartHorizon + weeksIn(pp) - 1;
    endIdx = startIdx + N_mhe - 1;
    
    Xholt_plot = plot(date(startIdx:endIdx), Sim.Holt.(currentWeekField).sim_yhat(6,:), 'LineWidth', 1.5);
    Xholt_color = Xholt_plot.Color;
    fill([date(startIdx:endIdx) flip(date(startIdx:endIdx))],...
         [Sim.Holt.(currentWeekField).sim_yhat_upper(6,:) flip(Sim.Holt.(currentWeekField).sim_yhat_lower(6,:))],...
         Xholt_color, 'FaceAlpha', .3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    hold on
    Xcost_plt = plot(date(startIdx:endIdx), Sim.Const.(currentWeekField)(6,:), 'LineWidth', 2,'LineStyle','-.');

    xline(date(startIdx), ":", 'HandleVisibility', 'off');
    
    hold on 
end
lgd = legend('Real Data', 'MHE results', sprintf('%d days Future Forecast', N_mhe), 'Interpreter', 'Latex');
lgd.FontSize = 18;
xlim([date(1) date(end)])
title('E-Expired', 'Interpreter', 'Latex'); 


%% Plotting of results (Prophet)

% I 
figure(21)
plot(date(1:end), I_data, '--') % Plot the real data
hold on
plot(date(N_mhe:358+N_mhe-1), Xsim.dyn(2,:), 'LineWidth', 1.5) % Plot the new fitted model

for pp = 1:length(weeksIn)
    currentWeekField = sprintf('wk%d', pp-1);
    startIdx = N_mhe + StartHorizon + weeksIn(pp) - 1;
    endIdx = startIdx + N_mhe - 1;
    
    XProph_plot = plot(date(startIdx:endIdx), Sim.Proph.(currentWeekField).sim_yhat(2,:), 'LineWidth', 1.5);
    XProph_color = XProph_plot.Color;
    fill([date(startIdx:endIdx) flip(date(startIdx:endIdx))],...
         [Sim.Proph.(currentWeekField).sim_yhat_upper(2,:) flip(Sim.Proph.(currentWeekField).sim_yhat_lower(2,:))],...
         XProph_color, 'FaceAlpha', .3, 'EdgeColor', 'none', 'HandleVisibility', 'off');

    xline(date(startIdx), ":", 'HandleVisibility', 'off');
    
    hold on 
end
lgd = legend('Real Data', 'MHE results', sprintf('%d days Future Forecast', N_mhe), 'Interpreter', 'Latex');
lgd.FontSize = 18;
xlim([date(1) date(end)])
title('I-Infected', 'Interpreter', 'Latex');


% D 
figure(22)
plot(date(1:end), D_data, '--') % Plot the real data
hold on
plot(date(N_mhe:358+N_mhe-1), Xsim.dyn(3,:), 'LineWidth', 1.5) % Plot the new fitted model

for pp = 1:length(weeksIn)
    currentWeekField = sprintf('wk%d', pp-1);
    startIdx = N_mhe + StartHorizon + weeksIn(pp) - 1;
    endIdx = startIdx + N_mhe - 1;
    
    XProph_plot = plot(date(startIdx:endIdx), Sim.Proph.(currentWeekField).sim_yhat(3,:), 'LineWidth', 1.5);
    XProph_color = XProph_plot.Color;
    fill([date(startIdx:endIdx) flip(date(startIdx:endIdx))],...
         [Sim.Proph.(currentWeekField).sim_yhat_upper(3,:) flip(Sim.Proph.(currentWeekField).sim_yhat_lower(3,:))],...
         XProph_color, 'FaceAlpha', .3, 'EdgeColor', 'none', 'HandleVisibility', 'off');

    xline(date(startIdx), ":", 'HandleVisibility', 'off');
    
    hold on
end
lgd = legend('Real Data', 'MHE results', sprintf('%d days Future Forecast', N_mhe), 'Interpreter', 'Latex');
lgd.FontSize = 18;
xlim([date(1) date(end)])
title('D-Infected', 'Interpreter', 'Latex'); 


% T 
figure(23)
plot(date(1:end), T_data, '--') % Plot the real data
hold on
plot(date(N_mhe:358+N_mhe-1), Xsim.dyn(4,:), 'LineWidth', 1.5) % Plot the new fitted model

for pp = 1:length(weeksIn)
    currentWeekField = sprintf('wk%d', pp-1);
    startIdx = N_mhe + StartHorizon + weeksIn(pp) - 1;
    endIdx = startIdx + N_mhe - 1;
    
    XProph_plot = plot(date(startIdx:endIdx), Sim.Proph.(currentWeekField).sim_yhat(4,:), 'LineWidth', 1.5);
    XProph_color = XProph_plot.Color;
    fill([date(startIdx:endIdx) flip(date(startIdx:endIdx))],...
         [Sim.Proph.(currentWeekField).sim_yhat_upper(4,:) flip(Sim.Proph.(currentWeekField).sim_yhat_lower(4,:))],...
         XProph_color, 'FaceAlpha', .3, 'EdgeColor', 'none', 'HandleVisibility', 'off');

    xline(date(startIdx), ":", 'HandleVisibility', 'off');
    
    hold on
end
lgd = legend('Real Data', 'MHE results', sprintf('%d days Future Forecast', N_mhe), 'Interpreter', 'Latex');
lgd.FontSize = 18;
xlim([date(1) date(end)])
title('T-Threatned', 'Interpreter', 'Latex'); 


% H 
figure(24)
plot(date(1:end), H_data, '--') % Plot the real data
hold on
plot(date(N_mhe:358+N_mhe-1), Xsim.dyn(5,:), 'LineWidth', 1.5) % Plot the new fitted model

for pp = 1:length(weeksIn)
    currentWeekField = sprintf('wk%d', pp-1);
    startIdx = N_mhe + StartHorizon + weeksIn(pp) - 1;
    endIdx = startIdx + N_mhe - 1;
    
    XProph_plot = plot(date(startIdx:endIdx), Sim.Proph.(currentWeekField).sim_yhat(5,:), 'LineWidth', 1.5);
    XProph_color = XProph_plot.Color;
    fill([date(startIdx:endIdx) flip(date(startIdx:endIdx))],...
         [Sim.Proph.(currentWeekField).sim_yhat_upper(5,:) flip(Sim.Proph.(currentWeekField).sim_yhat_lower(5,:))],...
         XProph_color, 'FaceAlpha', .3, 'EdgeColor', 'none', 'HandleVisibility', 'off');

    xline(date(startIdx), ":", 'HandleVisibility', 'off');
    
    hold on 
end
lgd = legend('Real Data', 'MHE results', sprintf('%d days Future Forecast', N_mhe), 'Interpreter', 'Latex');
lgd.FontSize = 18;
xlim([date(1) date(end)])
title('H-Healed', 'Interpreter', 'Latex'); 

% E - first week
figure(25)
plot(date(1:end), E_data, '--') % Plot the real data
hold on
plot(date(N_mhe:358+N_mhe-1), Xsim.dyn(6,:), 'LineWidth', 1.5) % Plot the new fitted model

for pp = 1:length(weeksIn)
    currentWeekField = sprintf('wk%d', pp-1);
    startIdx = N_mhe + StartHorizon + weeksIn(pp) - 1;
    endIdx = startIdx + N_mhe - 1;
    
    XProph_plot = plot(date(startIdx:endIdx), Sim.Proph.(currentWeekField).sim_yhat(6,:), 'LineWidth', 1.5);
    XProph_color = XProph_plot.Color;
    fill([date(startIdx:endIdx) flip(date(startIdx:endIdx))],...
         [Sim.Proph.(currentWeekField).sim_yhat_upper(6,:) flip(Sim.Proph.(currentWeekField).sim_yhat_lower(6,:))],...
         XProph_color, 'FaceAlpha', .3, 'EdgeColor', 'none', 'HandleVisibility', 'off');

    xline(date(startIdx), ":", 'HandleVisibility', 'off');
    
    hold on 
end
lgd = legend('Real Data', 'MHE results', sprintf('%d days Future Forecast', N_mhe), 'Interpreter', 'Latex');
lgd.FontSize = 18;
xlim([date(1) date(end)])
title('E-Expired', 'Interpreter', 'Latex'); 