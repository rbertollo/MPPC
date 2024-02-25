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

par = [0.6 0.3 0.02 0.03 0.05 0.09];        % Guessed coefficients to simulate the state dynamics
eqns_subs = subs(eqns, [alpha, gamma, delta, sigma, tau, lambda], par);

% Passing the vector of equations 
[~,~,f_ODE] = getdynamics(eqns_subs,vars,coefs);

x0 = [S_data(1), I_data(1), D_data(1), T_data(1), H_data(1), E_data(1)]; % Replace with actual initial conditions

% Define the time span for the simulation
tspan = 1:399; % Replace with the actual start and end times

% Solve the ODE for the ststem dynamics 
opts = odeset('RelTol',1e-3,'AbsTol',1e-6);
[t, x] = ode45(f_ODE, tspan, x0,opts);

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

matrixGSA1 = 
martixGSA2 =

%% Create a color-coded plot of the matrix 
shiftedColormap = colormap('copper') * 0.6 + 0.8; % This makes the darks lighter
shiftedColormap(shiftedColormap > 1) = 1;
shiftedColormap(shiftedColormap < 0) = 0;
text1 = {'\textit{S}', '\textit{I}', '\textit{D}', '\textit{T}', '\textit{H}', '\textit{E}'};
text2 = {'\textit{$\alpha$}', '\textit{$\gamma$}', '\textit{$\delta$}', '\textit{$\sigma$}', '\textit{$\tau$}', '\textit{$\lambda$}'};

figure(1)
imagesc(SensMat_sq)
colormap(shiftedColormap)
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