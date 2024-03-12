%% Python plots reported on matlab
clc
clear all
close all

load("/Users/marcodelloro/Desktop/MPPC/results.mat")
load('/Users/marcodelloro/Desktop/Thesis/MSc-Thesis-TUe/Data_Collection/Italian Dataset/SIDTTHE_data_DEF.mat');
% set(0,'DefaultFigureWindowStyle','docked');

dataset = readtable('/Users/marcodelloro/Desktop/Thesis/MSc-Thesis-TUe/Data_Collection/Italian Dataset/Italy_complete_dataset.xlsx');
load("results2.mat")
% load("resultsAUS.mat")

Npop = 59240329; % Total Population of Italy
N = 399;       % MPC horizon length
N_mhe = 21;

ForeWk0 = readtable('/Users/marcodelloro/Desktop/Thesis/MSc-Thesis-TUe/Forecasting/TrialPython/Forecasts21d_CI_test.xlsx', 'Sheet', 'Horizon_0_days');
ForeWk1 = readtable('/Users/marcodelloro/Desktop/Thesis/MSc-Thesis-TUe/Forecasting/TrialPython/Forecasts21d_CI_test.xlsx', 'Sheet', 'Horizon_7_days');
ForeWk2 = readtable('/Users/marcodelloro/Desktop/Thesis/MSc-Thesis-TUe/Forecasting/TrialPython/Forecasts21d_CI_test.xlsx', 'Sheet', 'Horizon_14_days');
ForeWk3 = readtable('/Users/marcodelloro/Desktop/Thesis/MSc-Thesis-TUe/Forecasting/TrialPython/Forecasts21d_CI_test.xlsx', 'Sheet', 'Horizon_21_days');
ForeWk4 = readtable('/Users/marcodelloro/Desktop/Thesis/MSc-Thesis-TUe/Forecasting/TrialPython/Forecasts21d_CI_test.xlsx', 'Sheet', 'Horizon_28_days');
ForeWk5 = readtable('/Users/marcodelloro/Desktop/Thesis/MSc-Thesis-TUe/Forecasting/TrialPython/Forecasts21d_CI_test.xlsx', 'Sheet', 'Horizon_35_days');


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

date = SIDTTHE_data{1,1}.date;

%% Results Plot
% Figure of the alpha trend related to policy In italy

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
% policy = {'Minimum restrictions', ''}
% handle = {'Minimum restrictions', ''}

for ii = 1:length(policy_dates)-1

    area.x(ii, :) = [policy_dates(ii) policy_dates(ii) policy_dates(ii+1) policy_dates(ii+1)];
    area.y_alpha(ii, :) = [0 max(results.par.alpha)*1.5 max(results.par.alpha)*1.5 0];
end
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
title('\textbf{$\alpha$ }','Interpreter','latex')
grid on
legend('Interpreter','latex','location','southeast')
xlim([date(1+N_mhe), date(end-N_mhe)])
ylim([0,max(results.par.alpha)*1.5])
set(gca, 'TickLabelInterpreter', 'Latex')

figure(6)
plot(date(N_mhe:N-N_mhe),results.par.gamma, 'LineWidth', 1.5, 'Color','k')
xlim([date(1+N_mhe), date(end-N_mhe)])
ylim([0,max(results.par.gamma)*1.5])
title('\textbf{$\gamma$}','Interpreter','latex')

figure(7)
plot(date(N_mhe:N-N_mhe),results.par.delta, 'LineWidth', 1.5, 'Color','k')
xlim([date(1+N_mhe), date(end-N_mhe)])
ylim([0,max(results.par.delta)*1.5])
title('\textbf{$\delta$}','Interpreter','latex')

figure(8)
plot(date(N_mhe:N-N_mhe),results.par.sigma, 'LineWidth', 1.5, 'Color','k')
xlim([date(1+N_mhe), date(end-N_mhe)])
ylim([0,max(results.par.sigma)*1.5])
title('\textbf{$\sigma$}','Interpreter','latex')

figure(9)
plot(date(N_mhe:N-N_mhe),results.par.tau, 'LineWidth', 1.5, 'Color','k')
xlim([date(1+N_mhe), date(end-N_mhe)])
ylim([0,max(results.par.tau)*1.5])
title('\textbf{$\tau$}','Interpreter','latex') 

%% Prediction Plot HOLT
dd = 35;

% Plotting of alpha - Holt
for pp = 1:length(weeksIn)
    currentWeekField = sprintf('wk%d', pp-1);
    startIdx = N_mhe + StartHorizon + weeksIn(pp) ;
    endIdx = startIdx + N_mhe;
    figure;
    
    plot(date(startIdx-dd:startIdx), results.par.alpha(startIdx-dd-N_mhe:startIdx-N_mhe),'LineWidth',2, 'DisplayName', '$\alpha$ Past')
    hold on
    Xholt_plot = plot(date(startIdx:endIdx), HoltFore.(currentWeekField).alpha_yhat,'LineWidth',2,'DisplayName', '$\alpha$ Forecasted');
    hold on
    plot(date(startIdx:endIdx), results.par.alpha(startIdx-N_mhe:endIdx-N_mhe),'k','LineWidth',2,'LineStyle','--','DisplayName', '$\alpha$ Future')
    hold on
    Xholt_color = Xholt_plot.Color;
    fill([date(startIdx:endIdx) flip(date(startIdx:endIdx))],...
         [HoltFore.(currentWeekField).alpha_yhat_upper' flip(HoltFore.(currentWeekField).alpha_yhat_lower)'],...
          Xholt_color, 'FaceAlpha', .3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    ylabel('Coefficients Values','Interpreter','latex')
    % title('\textbf{$\alpha$ }','Interpreter','latex')
    grid on
    % lgd = legend('Interpreter','latex','location','southeast');
    ylim([0,max(results.par.alpha)*1.5])
    set(gca, 'TickLabelInterpreter', 'Latex')
    box on 

end


% Plotting of Gamma - Holt
for pp = 1:length(weeksIn)
    currentWeekField = sprintf('wk%d', pp-1);
    startIdx = N_mhe + StartHorizon + weeksIn(pp) ;
    endIdx = startIdx + N_mhe;
    figure;
    
    plot(date(startIdx-dd:startIdx), results.par.gamma(startIdx-dd-N_mhe:startIdx-N_mhe),'LineWidth',1.5, 'DisplayName', '$\gamma$ Past')
    hold on
    Xholt_plot = plot(date(startIdx:endIdx), HoltFore.(currentWeekField).gamma_yhat,'LineWidth',1.5,'DisplayName', '$\gamma$ Forecasted');
    hold on
    plot(date(startIdx:endIdx), results.par.gamma(startIdx-N_mhe:endIdx-N_mhe),'k','LineWidth',1.5,'LineStyle','--','DisplayName', '$\gamma$ Future')
    hold on
    Xholt_color = Xholt_plot.Color;
    fill([date(startIdx:endIdx) flip(date(startIdx:endIdx))],...
         [HoltFore.(currentWeekField).gamma_yhat_upper' flip(HoltFore.(currentWeekField).gamma_yhat_lower)'],...
          Xholt_color, 'FaceAlpha', .3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    ylabel('Coefficients Values','Interpreter','latex')
    title('\textbf{$\gamma$ }','Interpreter','latex')
    grid on
    lgd = legend('Interpreter','latex','location','southeast');
    lgd.FontSize = 18;
    ylim([0,max(results.par.gamma)*1.5])
    set(gca, 'TickLabelInterpreter', 'Latex')

end


% Plotting of delta - Holt
for pp = 1:length(weeksIn)
    currentWeekField = sprintf('wk%d', pp-1);
    startIdx = N_mhe + StartHorizon + weeksIn(pp) ;
    endIdx = startIdx + N_mhe;
    figure;
    
    plot(date(startIdx-dd:startIdx), results.par.delta(startIdx-dd-N_mhe:startIdx-N_mhe),'LineWidth',2, 'DisplayName', '$\delta$ Past')
    hold on
    Xholt_plot = plot(date(startIdx:endIdx), HoltFore.(currentWeekField).delta_yhat,'LineWidth',2,'DisplayName', '$\delta$ Forecasted');
    hold on
    plot(date(startIdx:endIdx), results.par.delta(startIdx-N_mhe:endIdx-N_mhe),'k','LineWidth',2,'LineStyle','--','DisplayName', '$\delta$ Future')
    hold on
    Xholt_color = Xholt_plot.Color;
    fill([date(startIdx:endIdx) flip(date(startIdx:endIdx))],...
         [HoltFore.(currentWeekField).delta_yhat_upper' flip(HoltFore.(currentWeekField).delta_yhat_lower)'],...
          Xholt_color, 'FaceAlpha', .3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    ylabel('Coefficients Values','Interpreter','latex')
    % title('\textbf{$\delta$ }','Interpreter','latex')
    grid on
    % lgd = legend('Interpreter','latex','location','northeast');
    % lgd.FontSize = 18;
    set(gca, 'TickLabelInterpreter', 'Latex')
    box on 

end


% Plotting of sigma - Holt
for pp = 1:length(weeksIn)
    currentWeekField = sprintf('wk%d', pp-1);
    startIdx = N_mhe + StartHorizon + weeksIn(pp) ;
    endIdx = startIdx + N_mhe;
    figure;
    
    plot(date(startIdx-dd:startIdx), results.par.sigma(startIdx-dd-N_mhe:startIdx-N_mhe),'LineWidth',1.5, 'DisplayName', '$\sigma$ Past')
    hold on
    Xholt_plot = plot(date(startIdx:endIdx), HoltFore.(currentWeekField).sigma_yhat,'LineWidth',1.5,'DisplayName', '$\sigma$ Forecasted');
    hold on
    plot(date(startIdx:endIdx), results.par.sigma(startIdx-N_mhe:endIdx-N_mhe),'k','LineWidth',1.5,'LineStyle','--','DisplayName', '$\sigma$ Future')
    hold on
    Xholt_color = Xholt_plot.Color;
    fill([date(startIdx:endIdx) flip(date(startIdx:endIdx))],...
         [HoltFore.(currentWeekField).sigma_yhat_upper' flip(HoltFore.(currentWeekField).sigma_yhat_lower)'],...
          Xholt_color, 'FaceAlpha', .3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    ylabel('Coefficients Values','Interpreter','latex')
    title('\textbf{$\sigma$ }','Interpreter','latex')
    grid on
    lgd = legend('Interpreter','latex','location','southeast');
    lgd.FontSize = 18;
    ylim([0,max(results.par.sigma)*1.5])
    set(gca, 'TickLabelInterpreter', 'Latex')

end



% Plotting of tau - Holt
for pp = 1:length(weeksIn)
    currentWeekField = sprintf('wk%d', pp-1);
    startIdx = N_mhe + StartHorizon + weeksIn(pp) ;
    endIdx = startIdx + N_mhe;
    figure;
    
    plot(date(startIdx-dd:startIdx), results.par.tau(startIdx-dd-N_mhe:startIdx-N_mhe),'LineWidth',1.5, 'DisplayName', '$\tau$ Past')
    hold on
    Xholt_plot = plot(date(startIdx:endIdx), HoltFore.(currentWeekField).tau_yhat,'LineWidth',1.5,'DisplayName', '$\tau$ Forecasted');
    hold on
    plot(date(startIdx:endIdx), results.par.tau(startIdx-N_mhe:endIdx-N_mhe),'k','LineWidth',1.5,'LineStyle','--','DisplayName', '$\tau$ Future')
    hold on
    Xholt_color = Xholt_plot.Color;
    fill([date(startIdx:endIdx) flip(date(startIdx:endIdx))],...
         [HoltFore.(currentWeekField).tau_yhat_upper' flip(HoltFore.(currentWeekField).tau_yhat_lower)'],...
          Xholt_color, 'FaceAlpha', .3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    ylabel('Coefficients Values','Interpreter','latex')
    title('\textbf{$\alpha$ }','Interpreter','latex')
    grid on
    lgd = legend('Interpreter','latex','location','southeast');
    lgd.FontSize = 18;
    ylim([0,max(results.par.tau)*1.5])
    set(gca, 'TickLabelInterpreter', 'Latex')

end


%% Prediction Plot PROPHET
dd = 35;

% Plotting of alpha
for pp = 1:length(weeksIn)
    currentWeekField = sprintf('wk%d', pp-1);
    startIdx = N_mhe + StartHorizon + weeksIn(pp) ;
    endIdx = startIdx + N_mhe;
    figure;
    
    plot(date(startIdx-dd:startIdx), results.par.alpha(startIdx-dd-N_mhe:startIdx-N_mhe),'LineWidth',1.5, 'DisplayName', '$\alpha$ Past')
    hold on
    XProph_plot = plot(date(startIdx:endIdx), ProphFore.(currentWeekField).alpha_Prophet_yhat,'LineWidth',1.5,'DisplayName', '$\alpha$ Forecasted');
    hold on
    plot(date(startIdx:endIdx), results.par.alpha(startIdx-N_mhe:endIdx-N_mhe),'k','LineWidth',1.5,'LineStyle','--','DisplayName', '$\alpha$ Future')
    hold on
    XProph_color = XProph_plot.Color;
    fill([date(startIdx:endIdx) flip(date(startIdx:endIdx))],...
         [ProphFore.(currentWeekField).alpha_Prophet_yhat_upper' flip(ProphFore.(currentWeekField).alpha_Prophet_yhat_lower)'],...
          XProph_color, 'FaceAlpha', .3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    ylabel('Coefficients Values','Interpreter','latex')
    title('\textbf{$\alpha$ }','Interpreter','latex')
    grid on
    lgd = legend('Interpreter','latex','location','southeast');
    lgd.FontSize = 18;
    ylim([0,max(results.par.alpha)*1.5])
    set(gca, 'TickLabelInterpreter', 'Latex')

end


% Plotting of Gamma - Proph
for pp = 1:length(weeksIn)
    currentWeekField = sprintf('wk%d', pp-1);
    startIdx = N_mhe + StartHorizon + weeksIn(pp) ;
    endIdx = startIdx + N_mhe;
    figure;
    
    plot(date(startIdx-dd:startIdx), results.par.gamma(startIdx-dd-N_mhe:startIdx-N_mhe),'LineWidth',1.5, 'DisplayName', '$\gamma$ Past')
    hold on
    XProph_plot = plot(date(startIdx:endIdx), ProphFore.(currentWeekField).gamma_Prophet_yhat,'LineWidth',1.5,'DisplayName', '$\gamma$ Forecasted');
    hold on
    plot(date(startIdx:endIdx), results.par.gamma(startIdx-N_mhe:endIdx-N_mhe),'k','LineWidth',1.5,'LineStyle','--','DisplayName', '$\gamma$ Future')
    hold on
    XProph_color = XProph_plot.Color;
    fill([date(startIdx:endIdx) flip(date(startIdx:endIdx))],...
         [ProphFore.(currentWeekField).gamma_Prophet_yhat_upper' flip(ProphFore.(currentWeekField).gamma_Prophet_yhat_lower)'],...
          XProph_color, 'FaceAlpha', .3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    ylabel('Coefficients Values','Interpreter','latex')
    title('\textbf{$\gamma$ }','Interpreter','latex')
    grid on
    lgd = legend('Interpreter','latex','location','southeast');
    lgd.FontSize = 18;
    ylim([0,max(results.par.gamma)*1.5])
    set(gca, 'TickLabelInterpreter', 'Latex')

end


% Plotting of delta - Proph
for pp = 1:length(weeksIn)
    currentWeekField = sprintf('wk%d', pp-1);
    startIdx = N_mhe + StartHorizon + weeksIn(pp) ;
    endIdx = startIdx + N_mhe;
    figure;
    
    plot(date(startIdx-dd:startIdx), results.par.delta(startIdx-dd-N_mhe:startIdx-N_mhe),'LineWidth',1.5, 'DisplayName', '$\delta$ Past')
    hold on
    XProph_plot = plot(date(startIdx:endIdx), ProphFore.(currentWeekField).delta_Prophet_yhat,'LineWidth',1.5,'DisplayName', '$\delta$ Forecasted');
    hold on
    plot(date(startIdx:endIdx), results.par.delta(startIdx-N_mhe:endIdx-N_mhe),'k','LineWidth',1.5,'LineStyle','--','DisplayName', '$\delta$ Future')
    hold on
    XProph_color = XProph_plot.Color;
    fill([date(startIdx:endIdx) flip(date(startIdx:endIdx))],...
         [ProphFore.(currentWeekField).delta_Prophet_yhat_upper' flip(ProphFore.(currentWeekField).delta_Prophet_yhat_lower)'],...
          XProph_color, 'FaceAlpha', .3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    ylabel('Coefficients Values','Interpreter','latex')
    title('\textbf{$\delta$ }','Interpreter','latex')
    grid on
    lgd = legend('Interpreter','latex','location','southeast');
    lgd.FontSize = 18;
    set(gca, 'TickLabelInterpreter', 'Latex')

end


% Plotting of sigma - Proph
for pp = 1:length(weeksIn)
    currentWeekField = sprintf('wk%d', pp-1);
    startIdx = N_mhe + StartHorizon + weeksIn(pp) ;
    endIdx = startIdx + N_mhe;
    figure;
    
    plot(date(startIdx-dd:startIdx), results.par.sigma(startIdx-dd-N_mhe:startIdx-N_mhe),'LineWidth',1.5, 'DisplayName', '$\sigma$ Past')
    hold on
    XProph_plot = plot(date(startIdx:endIdx), ProphFore.(currentWeekField).sigma_Prophet_yhat,'LineWidth',1.5,'DisplayName', '$\sigma$ Forecasted');
    hold on
    plot(date(startIdx:endIdx), results.par.sigma(startIdx-N_mhe:endIdx-N_mhe),'k','LineWidth',1.5,'LineStyle','--','DisplayName', '$\sigma$ Future')
    hold on
    XProph_color = XProph_plot.Color;
    fill([date(startIdx:endIdx) flip(date(startIdx:endIdx))],...
         [ProphFore.(currentWeekField).sigma_Prophet_yhat_upper' flip(ProphFore.(currentWeekField).sigma_Prophet_yhat_lower)'],...
          XProph_color, 'FaceAlpha', .3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    ylabel('Coefficients Values','Interpreter','latex')
    title('\textbf{$\sigma$ }','Interpreter','latex')
    grid on
    lgd = legend('Interpreter','latex','location','southeast');
    lgd.FontSize = 18;
    ylim([0,max(results.par.sigma)*1.5])
    set(gca, 'TickLabelInterpreter', 'Latex')

end



% Plotting of tau - Proph
for pp = 1:length(weeksIn)
    currentWeekField = sprintf('wk%d', pp-1);
    startIdx = N_mhe + StartHorizon + weeksIn(pp) ;
    endIdx = startIdx + N_mhe;
    figure;
    
    plot(date(startIdx-dd:startIdx), results.par.tau(startIdx-dd-N_mhe:startIdx-N_mhe),'LineWidth',1.5, 'DisplayName', '$\tau$ Past')
    hold on
    XProph_plot = plot(date(startIdx:endIdx), ProphFore.(currentWeekField).tau_Prophet_yhat,'LineWidth',1.5,'DisplayName', '$\tau$ Forecasted');
    hold on
    plot(date(startIdx:endIdx), results.par.tau(startIdx-N_mhe:endIdx-N_mhe),'k','LineWidth',1.5,'LineStyle','--','DisplayName', '$\tau$ Future')
    hold on
    XProph_color = XProph_plot.Color;
    fill([date(startIdx:endIdx) flip(date(startIdx:endIdx))],...
         [ProphFore.(currentWeekField).tau_Prophet_yhat_upper' flip(ProphFore.(currentWeekField).tau_Prophet_yhat_lower)'],...
          XProph_color, 'FaceAlpha', .3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    ylabel('Coefficients Values','Interpreter','latex')
    title('\textbf{$\alpha$ }','Interpreter','latex')
    grid on
    lgd = legend('Interpreter','latex','location','southeast');
    lgd.FontSize = 18;
    ylim([0,max(results.par.tau)*1.5])
    set(gca, 'TickLabelInterpreter', 'Latex')

end
