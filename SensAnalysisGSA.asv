%% Sensitivity analysis on the SIDTHE model - GSA method

clear all
close all
clc
load("InitialCond.mat")

%% Generation of the LHS data

nSamples = 100;
minVals = [0.1 0.05 1e-4 1e-4 1e-4 0.05]; % Lower bound limit vector [alpha-gamma-delta-sigma-tau-lambda]
maxVals = [0.7 0.5 0.09 0.09 0.09 0.2]; % higher bound limit vector [alpha-gamma-delta-sigma-tau-lambda]

sampledValues3D = zeros(nSamples, length(minVals), 2);

for i = 1:2
    lhs = lhsdesign(nSamples, length(minVals)); 
    scaledLHS = minVals + lhs .* (maxVals - minVals);
    mat3D(:,:,i) = scaledLHS;
end

A = mat3D(:,:,1);
B = mat3D(:,:,2);

% Initialize struct to hold the combined matrices
C = struct();

for j = 1:6
    c_idx = sprintf('c%d', j);
    C.(c_idx) = A;
    C.(c_idx)(:, j) = B(:, j);
end

%% Simulation of the system with 100 different base parameters x2 times

syms S I D T H E % Variables
syms alpha gamma delta lambda sigma tau % Coefs

eqns = [    -S * (alpha * I) + lambda * I;...
             S * (alpha * I) - (gamma+lambda) * I;...
             I * gamma - D * (lambda + delta);...
             delta * D - (tau + sigma)*T;...
             lambda * (D) + T * sigma;...
             tau*T];

vars = [S I D T H E];
coefs = [alpha gamma delta sigma tau lambda];

y_A = zeros(399, 6, 100); % output storage matrix for A
y_B = zeros(399, 6, 100); % output storage matrix for B


for k = 1:2
    for ii = 1:size(A,1)
           par = mat3D(ii,:,k);
           eqns_subs = subs(eqns, [alpha, gamma, delta, sigma, tau, lambda], par);
           
           [~,~,f_ODE] = getdynamics(eqns_subs,vars,coefs);
           tspan = 1:399;
           opts = odeset('RelTol',1e-3,'AbsTol',1e-6);
           [t, x] = ode45(f_ODE, tspan, x0,opts);
           
           % Results storage 
           if k == 1
                y_A(:,:,ii) = x; 
           else
                y_B(:,:,ii) = x; 
           end
    end
end


% Computation of the Cj outputs Y_Cj
for j = 1:6
    c_idx = sprintf('c%d', j);
    for ii = 1:size(A,1)
        par = C.(c_idx)(ii, :);
        eqns_subs = subs(eqns, [alpha, gamma, delta, sigma, tau, lambda], par);
    
        [~,~,f_ODE] = getdynamics(eqns_subs, vars, coefs);
        tspan = 1:399;
        opts = odeset('RelTol',1e-3, 'AbsTol',1e-6);
        [t, x] = ode45(f_ODE, tspan, x0, opts);
        
        y_C.(c_idx)(:,:,ii) = x;
    end
end

% Calculate the Mean and Variance iot calculate sensitivity index (main effect, Si and total effect, STi)
Emain = mean(y_B, 3);
Etot = mean(y_A, 3).^2;

Vmain = mean((y_B).^2, 3) - Emain.^2;
Vtot = mean((y_A).^2, 3) - Etot.^2;

N=100;

% Compute the sensitivity index S_mainline for each element
for jj = 1:6
    c_idx = sprintf('c%d', jj);
    sumSqDiff = sum((y_B - y_C.(c_idx)).^2, 3);
    sumSqDiff_tot = sum((y_A - y_C.(c_idx)).^2, 3);

    Di = Vmain - sumSqDiff/(2*N);
    Di_tot = Vtot - sumSqDiff_tot/(2*N);

    s(:,jj,:) = Di'./Vmain';
    s_tot(:,jj,:) = Di_tot'./Vtot';
end

% Matlab code final results in L2 norm
s_l2 = sqrt(sum(s.^2,3)); % main sensitivity
stot_l2 = sqrt(sum(s_tot.^2,3)); % total sensitivity

% Doma/ Julia results 
s_l2doma = [ 0.495329  0.0188935   1.54295e-6  -3.23261e-7  -2.68334e-7  -0.0903141;
             0.212264  0.713847   -7.94141e-9  -1.68204e-8   1.80572e-7   0.0238113;
             0.586775  0.26544     0.0744124    3.48983e-8   1.24476e-6   0.224611;
             0.248215  0.221971    0.0272932   -0.0110097    0.0247515    0.0543262;
             0.654672  0.305631    0.0495066    0.0174521   -0.00227397   0.0786392;
             0.220084  0.128949   -0.0673783   -0.0274055   -0.00433011   0.10705       ];

stot_l2doma = [0.757088  0.52451   5.34127e-12  4.08903e-13  4.57446e-12  0.0722807;
           0.513223  1.00684   2.12281e-11  4.27393e-13  1.07463e-11  0.045626;
           0.846135  0.456466  0.0644594    2.47302e-12  6.59351e-11  0.296517;
           0.451054  0.353557  0.296036     0.146859     0.0707392    0.130103;
           0.807944  0.523788  0.0190629    0.00627838   0.00161554   0.054199;
           0.401785  0.251822  0.379464     0.0350034    0.117665     0.167419];

%% Create a color-coded plot of the matrix 
shiftedColormap = colormap('gray') * 0.5 + 0.65; % This makes the darks lighter
shiftedColormap(shiftedColormap > 1) = 1;
shiftedColormap(shiftedColormap < 0) = 0;
text1 = {'\textit{S}', '\textit{I}', '\textit{D}', '\textit{T}', '\textit{H}', '\textit{E}'};
text2 = {'\textit{$\alpha$}', '\textit{$\gamma$}', '\textit{$\delta$}', '\textit{$\sigma$}', '\textit{$\tau$}', '\textit{$\lambda$}'};


% Main sensitivity 
figure(1)
subplot(1,2,1)
imagesc(s_l2doma)
colormap(shiftedColormap)
set(gca, 'GridLineStyle', '-', 'LineWidth', 1, 'GridColor', 'k');
set(gca, 'XGrid', 'on', 'YGrid', 'on');
set(gca, 'Layer', 'top');

[m, n] = size(s_l2doma);
for i = 1:m
    for j = 1:n
      text(j, i, num2str(s_l2doma(i,j), '%0.2e'), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle',Interpreter='latex', FontSize=12);
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
title("Julia Results")

subplot(1,2,2)
imagesc(s_l2)
colormap(shiftedColormap)
set(gca, 'GridLineStyle', '-', 'LineWidth', 1, 'GridColor', 'k');
set(gca, 'XGrid', 'on', 'YGrid', 'on');
set(gca, 'Layer', 'top');

[m, n] = size(s_l2);
for i = 1:m
    for j = 1:n
      text(j, i, num2str(s_l2(i,j), '%0.2e'), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle',Interpreter='latex', FontSize=12);
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
title("MATLAB Results")


% Total sensitivity 
figure(2)
subplot(1,2,1)
imagesc(stot_l2)
colormap(shiftedColormap)
set(gca, 'GridLineStyle', '-', 'LineWidth', 1, 'GridColor', 'k');
set(gca, 'XGrid', 'on', 'YGrid', 'on');
set(gca, 'Layer', 'top');

[m, n] = size(stot_l2);
for i = 1:m
    for j = 1:n
      text(j, i, num2str(stot_l2(i,j), '%0.2e'), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle',Interpreter='latex', FontSize=12);
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
title("MATLAB Results")


subplot(1,2,2)
imagesc(stot_l2doma)
colormap(shiftedColormap)
set(gca, 'GridLineStyle', '-', 'LineWidth', 1, 'GridColor', 'k');
set(gca, 'XGrid', 'on', 'YGrid', 'on');
set(gca, 'Layer', 'top');

[m, n] = size(stot_l2doma);
for i = 1:m
    for j = 1:n
      text(j, i, num2str(stot_l2doma(i,j), '%0.2e'), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle',Interpreter='latex', FontSize=12);
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
title("Julia Results")



