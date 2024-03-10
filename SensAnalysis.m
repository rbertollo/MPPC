%% Sensitivity analysis on the SIDTHE model - O.I. Krivorotko method

clear all
close all
clc

% Data Loading for initial conditions

load('/Users/marcodelloro/Desktop/Thesis/MSc-Thesis-TUe/Data_Collection/Italian Dataset/SIDTTHE_data_DEF.mat');
Npop = 59240329; 

I_data = SIDTTHE_data{1,1}.data / Npop;     
D_data = SIDTTHE_data{2,1}.data / Npop;     
T1_data = SIDTTHE_data{3,1}.data / Npop;
T2_data = SIDTTHE_data{4,1}.data / Npop;
H_dataAug = SIDTTHE_data{6,1}.data / Npop;
H_data = SIDTTHE_data{7,1}.data / Npop;
E_data = SIDTTHE_data{5,1}.data  / Npop;
S_data = ones(length(I_data),1)' - (I_data + D_data + T1_data + T2_data + H_data + E_data );

T_data = T1_data + T2_data;

x0 = [S_data(1), I_data(1), D_data(1), T_data(1), H_data(1), E_data(1)]; % Replace with actual initial conditions

save('InitialCond.mat', 'x0')

%% Simulation of the dynamics of the system 

syms S I D T H E % Variables
syms alpha gamma delta lambda sigma tau % Coefs

eqns = [    -S * (alpha * I);...
             S * (alpha * I) - (gamma+lambda) * I;...
             I * gamma - D * (lambda + delta);...
             delta * D - (tau + sigma)*T;...
             lambda * (D) + T * sigma + lambda * I ;...
             tau*T];

vars = [S I D T H E];
coefs = [alpha gamma delta sigma tau lambda];

par = [0.35 0.1 0.02 0.03 0.05 0.09];        % Guessed coefficients to simulate the state dynamics
eqns_subs = subs(eqns, [alpha, gamma, delta, sigma, tau, lambda], par);

% Passing the vector of equations 
[~,~,f_ODE] = getdynamics(eqns_subs,vars,coefs);

x0 = [S_data(1), I_data(1), D_data(1), T_data(1), H_data(1), E_data(1)]; % Replace with actual initial conditions

% Define the time span for the simulation
tspan = 1:399; % Replace with the actual start and end times

% Solve the ODE for the ststem dynamics 
opts = odeset('RelTol',1e-3,'AbsTol',1e-6);
[t, x] = ode45(f_ODE, tspan, x0, opts);

J_x = jacobian(eqns_subs, vars);
J_q = jacobian(eqns, coefs);

% plot(SIDTTHE_data{1,1}.date(1:100), x, LineWidth=1.5)
% xlim([SIDTTHE_data{1,1}.date(1) SIDTTHE_data{1,1}.date(100)])


%% Solving of the ODE every timestep and Calculate Sensitivity Matrix

s0 = zeros(length(coefs),1 );
SensMat = [];
SensMat_sq = zeros(length(vars),length(coefs));

t = [1,1+1];

for ii=1:399
    
    J_xi = eval(subs(J_x, vars, x(ii,:)));
    J_qi = eval(subs(J_q, vars, x(ii,:)));
    
    for jj=1:length(coefs)
        
        [TS, s_out] = ode45(@(t, s) J_xi*s + J_qi(:,jj), t, s0);
        mtrx(:,jj) = s_out(end, :).';
        mtrx_sq = mtrx.^2;

    end
    
    SensMat_sq = SensMat_sq + mtrx_sq;
    SensMat = [SensMat; mtrx];
    SensMat = sqrt(SensMat);
end

matixGSA1 =[    0.477448  0.286754    9.89622e-8  4.47675e-8    3.39091e-8  0.0296979;
                 0.248988  0.425229    1.2657e-7   8.98342e-9    1.08108e-9  0.0510235;
                 0.35377   0.154062    0.0108664   7.8629e-8     2.04262e-9  0.16781;
                 0.22048   0.0923599   0.0599971   0.0217318     0.0181968   0.0983554;
                 0.438152  0.320464    0.00106837  0.000293899   9.70172e-5  0.0170487;
                 0.221364  0.0870397   0.0525733   0.0128336     0.0419382   0.0890415      ];

matrixGSA2 = [   0.682382  0.487687  2.36001e-11  2.2653e-12   3.89922e-12  0.065629;
                 0.465497  0.677296  7.20481e-11  6.62965e-12  7.94682e-12  0.158029;
                 0.614964  0.341343  0.0468327    7.25182e-12  1.11662e-11  0.388098;
                 0.559344  0.298792  0.214925     0.127201     0.123467     0.314248;
                 0.658182  0.538539  0.00271057   0.00200482   0.000685634  0.0460702;
                 0.575464  0.295631  0.223008     0.0623813    0.170996     0.314403        ];

%% Create a color-coded plot of the matrix 
SensMat_sq2 = SensMat_sq;
indices = [1 3; 1 4; 1 5; 2 3; 2 4; 2 5; 3 4; 3 5];

% Substituting 1e-15 into the specified cells
for i = 1:size(indices, 1)
    SensMat_sq2(indices(i, 1), indices(i, 2)) = 1e-12;
end

logData = log10(SensMat_sq2 + eps);

shiftedColormap = colormap('gray')* 0.45+ 0.5; % This makes the darks lighter
shiftedColormap(shiftedColormap > 1) = 1;
shiftedColormap(shiftedColormap < 0) = 0;
text1 = {'\textit{S}', '\textit{I}', '\textit{D}', '\textit{T}', '\textit{H}', '\textit{E}'};
text2 = {'\textit{$\alpha$}', '\textit{$\gamma$}', '\textit{$\delta$}', '\textit{$\sigma$}', '\textit{$\tau$}', '\textit{$\lambda$}'};

figure(1)
imagesc(logData)
colormap(shiftedColormap)
% contrastEnhancedData = histeq(normalizedData);
% imagesc(contrastEnhancedData); % Re-plot with enhanced contrast
set(gca, 'GridLineStyle', '-', 'LineWidth', 1, 'GridColor', 'k');
set(gca, 'XGrid', 'on', 'YGrid', 'on');
set(gca, 'Layer', 'top');

[m, n] = size(SensMat_sq);
for i = 1:m
    for j = 1:n
      text(j, i, num2str(SensMat_sq(i,j), '%0.2e'), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle',Interpreter='latex', FontSize=12);
      text(0, i, text1(i),'HorizontalAlignment', 'right','VerticalAlignment', 'middle','FontSize', 14, Interpreter='latex');
    end
end

for j = 1:n
    text(j, 0, text2(j), ...
         'HorizontalAlignment', 'center', ...
         'VerticalAlignment', 'top', ...
         'FontSize', 14, Interpreter='latex');
end

set(gca, 'XTick', []);
set(gca, 'YTick', []);

%% Sobol Matrix S1

logData = log10(martixGSA1 + eps);

shiftedColormap = colormap('gray')* 0.3 + 0.6; % This makes the darks lighter
shiftedColormap(shiftedColormap > 1) = 1;
shiftedColormap(shiftedColormap < 0) = 0;
text1 = {'\textit{S}', '\textit{I}', '\textit{D}', '\textit{T}', '\textit{H}', '\textit{E}'};
text2 = {'\textit{$\alpha$}', '\textit{$\gamma$}', '\textit{$\delta$}', '\textit{$\sigma$}', '\textit{$\tau$}', '\textit{$\lambda$}'};

figure(2)
imagesc(logData)
colormap(shiftedColormap)
% contrastEnhancedData = histeq(logData);
% imagesc(contrastEnhancedData); % Re-plot with enhanced contrast
set(gca, 'GridLineStyle', '-', 'LineWidth', 1, 'GridColor', 'k');
set(gca, 'XGrid', 'on', 'YGrid', 'on');
set(gca, 'Layer', 'top');

[m, n] = size(martixGSA1);
for i = 1:m
    for j = 1:n
      text(j, i, num2str(martixGSA1(i,j), '%0.2e'), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle',Interpreter='latex', FontSize=12);
      text(0, i, text1(i),'HorizontalAlignment', 'right','VerticalAlignment', 'middle','FontSize', 14, Interpreter='latex');
    end
end

for j = 1:n
    text(j, 0, text2(j), ...
         'HorizontalAlignment', 'center', ...
         'VerticalAlignment', 'top', ...
         'FontSize', 14, Interpreter='latex');
end

set(gca, 'XTick', []);
set(gca, 'YTick', []);

% Sobol Matrix Stot

logData2 = log10(matrixGSA2+eps);

shiftedColormap = colormap('gray') * 0.3 + 0.6; % This makes the darks lighter
shiftedColormap(shiftedColormap > 1) = 1;
shiftedColormap(shiftedColormap < 0) = 0;
text1 = {'\textit{S}', '\textit{I}', '\textit{D}', '\textit{T}', '\textit{H}', '\textit{E}'};
text2 = {'\textit{$\alpha$}', '\textit{$\gamma$}', '\textit{$\delta$}', '\textit{$\sigma$}', '\textit{$\tau$}', '\textit{$\lambda$}'};

figure(3)
imagesc(logData2)
colormap(shiftedColormap)
% contrastEnhancedData = histeq(logData2);
% imagesc(contrastEnhancedData); % Re-plot with enhanced contrast
set(gca, 'GridLineStyle', '-', 'LineWidth', 1, 'GridColor', 'k');
set(gca, 'XGrid', 'on', 'YGrid', 'on');
set(gca, 'Layer', 'top');

[m, n] = size(matrixGSA2);
for i = 1:m
    for j = 1:n
      text(j, i, num2str(matrixGSA2(i,j), '%0.2e'), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle',Interpreter='latex', FontSize=12);
      text(0, i, text1(i),'HorizontalAlignment', 'right','VerticalAlignment', 'middle','FontSize', 14, Interpreter='latex');
    end
end

for j = 1:n
    text(j, 0, text2(j), ...
         'HorizontalAlignment', 'center', ...
         'VerticalAlignment', 'top', ...
         'FontSize', 14, Interpreter='latex');
end

set(gca, 'XTick', []);
set(gca, 'YTick', []);

%% Histogram Bar plots

paramNames = {'$\alpha$', '$\gamma$', '$\delta$', '$\sigma$', '$\tau$', '$\lambda$'};
statesNames = {'\textit{\textbf{S}}', '\textit{\textbf{I}}', '\textit{\textbf{D}}', '\textit{\textbf{T}}',...
                '\textit{\textbf{H}}', '\textit{\textbf{E}}'};

figure(10)
for i = 1:6
    subplot(3, 2, i);
    bar([matrixGSA1(i, :); matrixGSA2(i, :)]', 'grouped');
    title((statesNames(i)), 'Interpreter', 'latex');
    ylabel('Sobol Index', 'Interpreter', 'latex');
    
    % Set the x-axis ticks to the parameter names
    set(gca, 'XTick', 1:6, 'XTickLabel', paramNames, 'TickLabelInterpreter', 'latex');
    
    % Only show legend on the first subplot for clarity
    if i == 1
        legend({'1st Order Effect', 'Total Effect'}, 'Interpreter', 'latex', 'Location', 'best');
    end
end