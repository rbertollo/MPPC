%% Script with Forecast of 5 weeks - Last version
clc
clear all 
close all

load('/Users/marcodelloro/Desktop/Thesis/MSc-Thesis-TUe/Data_Collection/Italian Dataset/SIDTTHE_data_DEF.mat');
% set(0,'DefaultFigureWindowStyle','docked');

dataset = readtable('/Users/marcodelloro/Desktop/Thesis/MSc-Thesis-TUe/Data_Collection/Italian Dataset/Italy_complete_dataset.xlsx');
dataset = dataset(190:588,:);
dataset.data = datetime(dataset.data, 'InputFormat', 'dd-MMM-yyyy HH:mm:ss');
load('/Users/marcodelloro/Desktop/Thesis/MSc-Thesis-TUe/Data_Collection/Italian Dataset/Var_Infected.mat')
load("results2.mat")

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

ForeWk0 = readtable('/Users/marcodelloro/Desktop/Thesis/MSc-Thesis-TUe/Forecasting/TrialPython/Forecasts21d_CI3.xlsx', 'Sheet', 'Horizon_0_days');
ForeWk1 = readtable('/Users/marcodelloro/Desktop/Thesis/MSc-Thesis-TUe/Forecasting/TrialPython/Forecasts21d_CI3.xlsx', 'Sheet', 'Horizon_7_days');
ForeWk2 = readtable('/Users/marcodelloro/Desktop/Thesis/MSc-Thesis-TUe/Forecasting/TrialPython/Forecasts21d_CI3.xlsx', 'Sheet', 'Horizon_14_days');
ForeWk3 = readtable('/Users/marcodelloro/Desktop/Thesis/MSc-Thesis-TUe/Forecasting/TrialPython/Forecasts21d_CI3.xlsx', 'Sheet', 'Horizon_21_days');
ForeWk4 = readtable('/Users/marcodelloro/Desktop/Thesis/MSc-Thesis-TUe/Forecasting/TrialPython/Forecasts21d_CI3.xlsx', 'Sheet', 'Horizon_28_days');
ForeWk5 = readtable('/Users/marcodelloro/Desktop/Thesis/MSc-Thesis-TUe/Forecasting/TrialPython/Forecasts21d_CI3.xlsx', 'Sheet', 'Horizon_35_days');

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

StartHorizon = 10*N_mhe; % we are here
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
                   results.stszero.E(StartHorizon+j+weeksIn(pp))           ];
                
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
                   results.stszero.E(StartHorizon+j+weeksIn(pp))          ];
                
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
                   results.stszero.E(StartHorizon+j+weeksIn(pp))      ];
                
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
plot(date(1:end), I_data, '--','LineWidth', 1.5,Color=[0, 0.4470, 0.7410]) % Plot the real data
colorsct=[0, 0.4470, 0.7410];
hold on
fill(var_area.x,var_area.Ivar/Npop, colorsct, 'FaceAlpha', .1, 'EdgeColor', 'none');
hold on
plot(date(N_mhe:N_mhe+StartHorizon-1), Xsim.dyn(2,1:StartHorizon), 'LineWidth', 2,'Color',[0.8500 0.3250 0.0980]) % Plot the new fitted model
hold on
plot(date(N_mhe+StartHorizon:358+N_mhe-1), Xsim.dyn(2,StartHorizon+1:end), 'LineWidth', 2,'Color',[0.94 0.73 0.6392], 'HandleVisibility', 'off')
hold on 
% fill(var_area.x, var_area.y/Npop, [0.5, 0.5, 0.5],'FaceAlpha', .3, 'EdgeColor', 'none', 'HandleVisibility', 'off')
for pp = 1:length(weeksIn)
    currentWeekField = sprintf('wk%d', pp-1);
    startIdx = N_mhe + StartHorizon + weeksIn(pp) - 1;
    endIdx = startIdx + N_mhe - 1;
    
    Xholt_plot = plot(date(startIdx:endIdx), Sim.Holt.(currentWeekField).sim_yhat(2,:), 'LineWidth', 2);
    Xholt_color = Xholt_plot.Color;
    fill([date(startIdx:endIdx) flip(date(startIdx:endIdx))],...
         [Sim.Holt.(currentWeekField).sim_yhat_upper(2,:) flip(Sim.Holt.(currentWeekField).sim_yhat_lower(2,:))],...
         Xholt_color, 'FaceAlpha', .3, 'EdgeColor', 'none');
    hold on
    % Xcost_plt = plot(date(startIdx:endIdx), Sim.Const.(currentWeekField)(2,:), 'LineWidth', 2,'LineStyle','-.');
    xline(date(startIdx), ":",'LineWidth', 1.5,'HandleVisibility', 'off');
    hold on 
end
lgd = legend('Real Data', 'Data 95\% confidence interval', 'MHE Fitted Data', sprintf('%d days Ahead Forecast', N_mhe), 'Forecast 95\% confidence interval', 'Interpreter', 'Latex');
% lgd.FontSize = 18;
% title('I-Infected', 'Interpreter', 'Latex');
yax = ylabel('Normalized Population','Interpreter','latex');
xlim([date(1+N_mhe), date(end-N_mhe)])
grid on
set(gca, 'TickLabelInterpreter', 'Latex')
box on 


% D 
figure(12)
plot(date(1:end), D_data, '--','LineWidth', 1.5) % Plot the real data
hold on
plot(date(N_mhe:N_mhe+StartHorizon-1), Xsim.dyn(3,1:StartHorizon), 'LineWidth', 2,'Color',[0.8500 0.3250 0.0980]) % Plot the new fitted model
hold on
plot(date(N_mhe+StartHorizon:358+N_mhe-1), Xsim.dyn(3,StartHorizon+1:end), 'LineWidth', 2,'Color',[0.94 0.73 0.6392], 'HandleVisibility', 'off')
hold on 
for pp = 1:length(weeksIn)
    currentWeekField = sprintf('wk%d', pp-1);
    startIdx = N_mhe + StartHorizon + weeksIn(pp) - 1;
    endIdx = startIdx + N_mhe - 1;
    
    Xholt_plot = plot(date(startIdx:endIdx), Sim.Holt.(currentWeekField).sim_yhat(3,:), 'LineWidth', 2);
    Xholt_color = Xholt_plot.Color;
    fill([date(startIdx:endIdx) flip(date(startIdx:endIdx))],...
         [Sim.Holt.(currentWeekField).sim_yhat_upper(3,:) flip(Sim.Holt.(currentWeekField).sim_yhat_lower(3,:))],...
         Xholt_color, 'FaceAlpha', .3, 'EdgeColor', 'none');
    hold on
    % Xcost_plt = plot(date(startIdx:endIdx), Sim.Const.(currentWeekField)(3,:), 'LineWidth', 2,'LineStyle','-.');


    xline(date(startIdx), ":", 'HandleVisibility', 'off');
    
    hold on
end
lgd = legend('Real Data', 'MHE Fitted Data', sprintf('%d days Ahead Forecast', N_mhe), 'Forecast 95\% confidence interval', 'Interpreter', 'Latex',Location='southeast'); %'northeast'
% lgd.FontSize = 18;
yax = ylabel('Normalized Population','Interpreter','latex');
xlim([date(100), date(end-N_mhe-80)])
grid on
set(gca, 'TickLabelInterpreter', 'Latex')
box on 


% T 
figure(13)
plot(date(1:end), T_data, '--','LineWidth', 1.5) % Plot the real data
hold on
plot(date(N_mhe:N_mhe+StartHorizon-1), Xsim.dyn(4,1:StartHorizon), 'LineWidth', 2,'Color',[0.8500 0.3250 0.0980]) % Plot the new fitted model
hold on
plot(date(N_mhe+StartHorizon-1:358+N_mhe-1), Xsim.dyn(4,StartHorizon:end), 'LineWidth', 2,'Color',[0.94 0.73 0.6392], 'HandleVisibility', 'off')
hold on 
for pp = 1:length(weeksIn)
    currentWeekField = sprintf('wk%d', pp-1);
    startIdx = N_mhe + StartHorizon + weeksIn(pp) - 1;
    endIdx = startIdx + N_mhe - 1;
    
    Xholt_plot = plot(date(startIdx:endIdx), Sim.Holt.(currentWeekField).sim_yhat(4,:), 'LineWidth', 2);
    Xholt_color = Xholt_plot.Color;
    fill([date(startIdx:endIdx) flip(date(startIdx:endIdx))],...
         [Sim.Holt.(currentWeekField).sim_yhat_upper(4,:) flip(Sim.Holt.(currentWeekField).sim_yhat_lower(4,:))],...
         Xholt_color, 'FaceAlpha', .3, 'EdgeColor', 'none');
    hold on
    % Xcost_plt = plot(date(startIdx:endIdx), Sim.Const.(currentWeekField)(4,:), 'LineWidth', 2,'LineStyle','-.');

    xline(date(startIdx), ":", 'HandleVisibility', 'off');
    
    hold on
end
% lgd = legend('Real Data', 'MHE Fitted Data', sprintf('%d days Ahead Forecast', N_mhe), 'Forecast 95\% confidence interval', 'Interpreter', 'Latex',Location='southeast');
% lgd.FontSize = 18;
yax = ylabel('Normalized Population','Interpreter','latex');
xlim([date(100), date(end-N_mhe-80)])
grid on
set(gca, 'TickLabelInterpreter', 'Latex')
box on


% % H 
% figure(14)
% plot(date(1:end), H_data, '--') % Plot the real data
% hold on
% plot(date(N_mhe:358+N_mhe-1), Xsim.dyn(5,:), 'LineWidth', 1.5) % Plot the new fitted model
% 
% for pp = 1:length(weeksIn)
%     currentWeekField = sprintf('wk%d', pp-1);
%     startIdx = N_mhe + StartHorizon + weeksIn(pp) - 1;
%     endIdx = startIdx + N_mhe - 1;
% 
%     Xholt_plot = plot(date(startIdx:endIdx), Sim.Holt.(currentWeekField).sim_yhat(5,:), 'LineWidth', 1.5);
%     Xholt_color = Xholt_plot.Color;
%     fill([date(startIdx:endIdx) flip(date(startIdx:endIdx))],...
%          [Sim.Holt.(currentWeekField).sim_yhat_upper(5,:) flip(Sim.Holt.(currentWeekField).sim_yhat_lower(5,:))],...
%          Xholt_color, 'FaceAlpha', .3, 'EdgeColor', 'none');
%     hold on
%     Xcost_plt = plot(date(startIdx:endIdx), Sim.Const.(currentWeekField)(5,:), 'LineWidth', 2,'LineStyle','-.');
% 
%     xline(date(startIdx), ":", 'HandleVisibility', 'off');
% 
%     hold on 
% end
% lgd = legend('Real Data', 'MHE results', sprintf('%d days Future Forecast', N_mhe), 'Interpreter', 'Latex');
% lgd.FontSize = 18;
% xlim([date(70), date(end-N_mhe-60)])
% title('H-Healed', 'Interpreter', 'Latex'); 

% E - first week
figure(15)
plot(date(1:end), E_data, '--', LineWidth=1.5) % Plot the real data
hold on
plot(date(N_mhe:N_mhe+StartHorizon-1), Xsim.dyn(6,1:StartHorizon), 'LineWidth', 2,'Color',[0.8500 0.3250 0.0980]) % Plot the new fitted model
hold on
plot(date(N_mhe+StartHorizon-1:358+N_mhe-1), Xsim.dyn(6,StartHorizon:end), 'LineWidth', 2,'Color',[0.94 0.73 0.6392], 'HandleVisibility', 'off')
hold on 
for pp = 1:length(weeksIn)
    currentWeekField = sprintf('wk%d', pp-1);
    startIdx = N_mhe + StartHorizon + weeksIn(pp) - 1;
    endIdx = startIdx + N_mhe - 1;
    
    Xholt_plot = plot(date(startIdx:endIdx), Sim.Holt.(currentWeekField).sim_yhat(6,:), 'LineWidth', 2);
    Xholt_color = Xholt_plot.Color;
    fill([date(startIdx:endIdx) flip(date(startIdx:endIdx))],...
         [Sim.Holt.(currentWeekField).sim_yhat_upper(6,:) flip(Sim.Holt.(currentWeekField).sim_yhat_lower(6,:))],...
         Xholt_color, 'FaceAlpha', .3, 'EdgeColor', 'none');
    hold on
    % Xcost_plt = plot(date(startIdx:endIdx), Sim.Const.(currentWeekField)(6,:), 'LineWidth', 2,'LineStyle','-.');

    xline(date(startIdx), ":", 'HandleVisibility', 'off');
    
    hold on 
end
% lgd = legend('Real Data', 'MHE Fitted Data', sprintf('%d days Ahead Forecast', N_mhe), 'Forecast 95\% confidence interval', 'Interpreter', 'Latex',Location='southeast');
% lgd.FontSize = 18;
yax = ylabel('Normalized Population','Interpreter','latex');
xlim([date(100), date(end-N_mhe-80)])
grid on
set(gca, 'TickLabelInterpreter', 'Latex')
box on 


%% Plotting of results (Prophet)

% I 
figure(21)
plot(date(1:end), I_data, '--','LineWidth', 1.5) % Plot the real data
hold on
plot(date(N_mhe:N_mhe+StartHorizon-1), Xsim.dyn(2,1:StartHorizon), 'LineWidth', 2,'Color',[0.8500 0.3250 0.0980]) % Plot the new fitted model
hold on
plot(date(N_mhe+StartHorizon:358+N_mhe-1), Xsim.dyn(2,StartHorizon+1:end), 'LineWidth', 2,'Color',[0.94 0.73 0.6392], 'HandleVisibility', 'off')
hold on 
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
xlim([date(100), date(end-N_mhe-80)])
title('I-Infected', 'Interpreter', 'Latex');


% D 
figure(22)
plot(date(1:end), D_data, '--','LineWidth', 1.5) % Plot the real data
hold on
plot(date(N_mhe:N_mhe+StartHorizon-1), Xsim.dyn(3,1:StartHorizon), 'LineWidth', 2,'Color',[0.8500 0.3250 0.0980]) % Plot the new fitted model
hold on
plot(date(N_mhe+StartHorizon-1:358+N_mhe-1), Xsim.dyn(3,StartHorizon:end), 'LineWidth', 2,'Color',[0.94 0.73 0.6392], 'HandleVisibility', 'off')
hold on 
for pp = 1:length(weeksIn)
    currentWeekField = sprintf('wk%d', pp-1);
    startIdx = N_mhe + StartHorizon + weeksIn(pp) - 1;
    endIdx = startIdx + N_mhe - 1;
    
    XProph_plot = plot(date(startIdx:endIdx), Sim.Proph.(currentWeekField).sim_yhat(3,:), 'LineWidth', 1.5);
    XProph_color = XProph_plot.Color;
    fill([date(startIdx:endIdx) flip(date(startIdx:endIdx))],...
         [Sim.Proph.(currentWeekField).sim_yhat_upper(3,:) flip(Sim.Proph.(currentWeekField).sim_yhat_lower(3,:))],...
         XProph_color, 'FaceAlpha', .3, 'EdgeColor', 'none');

    xline(date(startIdx), ":", 'HandleVisibility', 'off');
    
    hold on
end
% lgd = legend('Real Data', 'MHE Fitted Data', sprintf('%d days Ahead Forecast', N_mhe), 'Forecast 95\% confidence interval', 'Interpreter', 'Latex',Location='northeast');
% lgd.FontSize = 18;
yax = ylabel('Normalized Population','Interpreter','latex');
xlim([date(100), date(end-N_mhe-80)])
% ylim([0 0.021])
grid on
set(gca, 'TickLabelInterpreter', 'Latex')
box on 



% T 
figure(23)
plot(date(1:end), T_data, '--','LineWidth', 1.5) % Plot the real data
hold on
plot(date(N_mhe:N_mhe+StartHorizon-1), Xsim.dyn(4,1:StartHorizon), 'LineWidth', 2,'Color',[0.8500 0.3250 0.0980]) % Plot the new fitted model
hold on
plot(date(N_mhe+StartHorizon:358+N_mhe-1), Xsim.dyn(4,StartHorizon+1:end), 'LineWidth', 2,'Color',[0.94 0.73 0.6392], 'HandleVisibility', 'off')
hold on 

for pp = 1:length(weeksIn)
    currentWeekField = sprintf('wk%d', pp-1);
    startIdx = N_mhe + StartHorizon + weeksIn(pp) - 1;
    endIdx = startIdx + N_mhe - 1;
    
    XProph_plot = plot(date(startIdx:endIdx), Sim.Proph.(currentWeekField).sim_yhat(4,:), 'LineWidth', 1.5);
    XProph_color = XProph_plot.Color;
    fill([date(startIdx:endIdx) flip(date(startIdx:endIdx))],...
         [Sim.Proph.(currentWeekField).sim_yhat_upper(4,:) flip(Sim.Proph.(currentWeekField).sim_yhat_lower(4,:))],...
         XProph_color, 'FaceAlpha', .3, 'EdgeColor', 'none');

    xline(date(startIdx), ":", 'HandleVisibility', 'off');
    
    hold on
end
% lgd = legend('Real Data', 'MHE Fitted Data', sprintf('%d days Ahead Forecast', N_mhe), 'Forecast 95\% confidence interval', 'Interpreter', 'Latex',Location='northeast');
% lgd.FontSize = 18;
yax = ylabel('Normalized Population','Interpreter','latex');
xlim([date(100), date(end-N_mhe-80)])
grid on
set(gca, 'TickLabelInterpreter', 'Latex')
box on 

% % H 
% figure(24)
% plot(date(1:end), H_data, '--', 1.5) % Plot the real data
% hold on
% plot(date(N_mhe:358+N_mhe-1), Xsim.dyn(5,:), 'LineWidth', 1.5) % Plot the new fitted model
% 
% for pp = 1:length(weeksIn)
%     currentWeekField = sprintf('wk%d', pp-1);
%     startIdx = N_mhe + StartHorizon + weeksIn(pp) - 1;
%     endIdx = startIdx + N_mhe - 1;
% 
%     XProph_plot = plot(date(startIdx:endIdx), Sim.Proph.(currentWeekField).sim_yhat(5,:), 'LineWidth', 1.5);
%     XProph_color = XProph_plot.Color;
%     fill([date(startIdx:endIdx) flip(date(startIdx:endIdx))],...
%          [Sim.Proph.(currentWeekField).sim_yhat_upper(5,:) flip(Sim.Proph.(currentWeekField).sim_yhat_lower(5,:))],...
%          XProph_color, 'FaceAlpha', .3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
% 
%     xline(date(startIdx), ":", 'HandleVisibility', 'off');
% 
%     hold on 
% end
% lgd = legend('Real Data', 'MHE results', sprintf('%d days Future Forecast', N_mhe), 'Interpreter', 'Latex');
% lgd.FontSize = 18;
% xlim([date(70), date(end-N_mhe-60)])
% title('H-Healed', 'Interpreter', 'Latex'); 

% E - first week
figure(25)
plot(date(1:end), E_data, '--','LineWidth', 1.5) % Plot the real data
hold on
plot(date(N_mhe:N_mhe+StartHorizon-1), Xsim.dyn(6,1:StartHorizon), 'LineWidth', 2,'Color',[0.8500 0.3250 0.0980]) % Plot the new fitted model
hold on
plot(date(N_mhe+StartHorizon:358+N_mhe-1), Xsim.dyn(6,StartHorizon+1:end), 'LineWidth', 2,'Color',[0.94 0.73 0.6392], 'HandleVisibility', 'off')
hold on 

for pp = 1:length(weeksIn)
    currentWeekField = sprintf('wk%d', pp-1);
    startIdx = N_mhe + StartHorizon + weeksIn(pp) - 1;
    endIdx = startIdx + N_mhe - 1;
    
    XProph_plot = plot(date(startIdx:endIdx), Sim.Proph.(currentWeekField).sim_yhat(6,:), 'LineWidth', 1.5);
    XProph_color = XProph_plot.Color;
    fill([date(startIdx:endIdx) flip(date(startIdx:endIdx))],...
         [Sim.Proph.(currentWeekField).sim_yhat_upper(6,:) flip(Sim.Proph.(currentWeekField).sim_yhat_lower(6,:))],...
         XProph_color, 'FaceAlpha', .3, 'EdgeColor', 'none');

    xline(date(startIdx), ":", 'HandleVisibility', 'off');
    
    hold on 
end
% lgd = legend('Real Data', 'MHE Fitted Data', sprintf('%d days Ahead Forecast', N_mhe), 'Forecast 95\% confidence interval', 'Interpreter', 'Latex',Location='southeast');
% lgd.FontSize = 18;
yax = ylabel('Normalized Population','Interpreter','latex');
xlim([date(100), date(end-N_mhe-80)])
grid on
set(gca, 'TickLabelInterpreter', 'Latex')
box on 

%% Score Function 

% Score on E-deceased population

for pp = 1:1
    startIdx = N_mhe + StartHorizon + weeksIn(pp) - 1;
    startIdxMHE = StartHorizon + weeksIn(pp);

    currentWeekField = sprintf('wk%d', pp-1);
    weekend = [7,14,21];
    
    for ii = 1:3
        % score week by week for holt model
        SE_holt = abs( ( Sim.Holt.(currentWeekField).sim_yhat(6,weekend(ii)) - E_data(startIdx + weekend(ii)) ) ./ ( E_data(startIdx + weekend(ii)) -  E_data(startIdx) ) ) ;
        SE_holtMHE = abs( ( Sim.Holt.(currentWeekField).sim_yhat(6,weekend(ii)) - Xsim.dyn(6,startIdxMHE) ) ./ ( Xsim.dyn(6,startIdxMHE + weekend(ii)) -  Xsim.dyn(6,startIdxMHE) ) );
        SE_ProphMHE = abs( ( Sim.Proph.(currentWeekField).sim_yhat(6,weekend(ii)) - Xsim.dyn(6,startIdxMHE + weekend(ii)) ) ./ ( Xsim.dyn(6,startIdxMHE + weekend(ii)) -  Xsim.dyn(6,startIdxMHE) ) );
        score.holtData.(currentWeekField).(['week', num2str(weekend(ii))]) = SE_holt;
        score.holtMHE.(currentWeekField).(['week', num2str(weekend(ii))]) = SE_holtMHE;
        score.ProphMHE.(currentWeekField).(['week', num2str(weekend(ii))]) = SE_ProphMHE;

        % baseline score
        SE_baseline = abs( ( 2*E_data(startIdx) - E_data(startIdx - weekend(1)) - E_data(startIdx + weekend(ii)) ) ./ ( E_data(startIdx + weekend(ii)) -  E_data(startIdx) ) );
        SE_BaseMHE = abs( ( 2*Xsim.dyn(6,startIdxMHE) - Xsim.dyn(6,startIdxMHE - weekend(1)) - Xsim.dyn(6,startIdxMHE + weekend(ii)) ) ./ ( Xsim.dyn(6,startIdxMHE + weekend(ii)) -  Xsim.dyn(6,startIdxMHE) ) );
        score.baseData.(currentWeekField).(['week', num2str(weekend(ii))]) = SE_baseline;
        score.baseMHE.(currentWeekField).(['week', num2str(weekend(ii))]) = SE_BaseMHE;
    end
end

%% MAPE + RMSE scores on the compartments

% Gather the weeks to investigate from both data Sim.dyn

startIdx2 = StartHorizon + weeksIn(1);
endIdx2 = StartHorizon + weeksIn(end) + N_mhe - 1;

eval_period = Xsim.dyn(:,startIdx2:endIdx2);

% Divide the evaluation period in weeks
wks_mat= cell(1, 8);
for i = 1:8
    start_col = (i - 1) * 7 + 1;
    end_col = i * 7;
    wks_mat{i} = eval_period(:, start_col:end_col);
end

foreWks = [1,7,14];

for pp = 1:length(weeksIn)
    currentWeekField = sprintf('wk%d', pp-1);

    if strcmp(currentWeekField, 'wk0')
            wk1 = wks_mat{1};
            wk2 = wks_mat{2};
            wk3 = wks_mat{3};
        elseif strcmp(currentWeekField, 'wk1')
            wk1 = wks_mat{2};
            wk2 = wks_mat{3};
            wk3 = wks_mat{4};
        elseif strcmp(currentWeekField, 'wk2')
            wk1 = wks_mat{3};
            wk2 = wks_mat{4};
            wk3 = wks_mat{5};
        elseif strcmp(currentWeekField, 'wk3')
            wk1 = wks_mat{4};
            wk2 = wks_mat{5};
            wk3 = wks_mat{6};
        elseif strcmp(currentWeekField, 'wk4')
            wk1 = wks_mat{5};
            wk2 = wks_mat{6};
            wk3 = wks_mat{7};
        elseif strcmp(currentWeekField, 'wk5')
            wk1 = wks_mat{6};
            wk2 = wks_mat{7};
            wk3 = wks_mat{8};
    end

    for jj = 1:length(foreWks)
       y_pred1 = Sim.Holt.(currentWeekField).sim_yhat(:,foreWks(jj):foreWks(jj)+6);
       y_pred2 = Sim.Proph.(currentWeekField).sim_yhat(:,foreWks(jj):foreWks(jj)+6);

        % select the week
        if jj == 1
            selectedWk = wk1;
        elseif jj == 2
            selectedWk = wk2;
        elseif jj == 3
            selectedWk = wk3;
        end

        WeeklyMAPE = mapeFunc(selectedWk', y_pred1');
        WeeklyMAPE2 = mapeFunc(selectedWk', y_pred2');
        predictionWeekField = sprintf('wk%d', jj);
        MAPE_ResultsHolt.(currentWeekField).(predictionWeekField) = WeeklyMAPE;
        MAPE_ResultsProph.(currentWeekField).(predictionWeekField) = WeeklyMAPE2;

    end
end

