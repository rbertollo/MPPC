clc; clear; close all;
default_paper;

%% Load data
dataset = readtable('Italy_complete_dataset.xlsx');
dataset = dataset(190:588,:);
dataset.data = datetime(dataset.data, 'InputFormat', 'dd-MMM-yyyy HH:mm:ss');
% policy_idx = [1 40 69 116 141 190 242 294 399];
Npop = 59240329;

I_uf = (dataset.ricoverati + ...
            dataset.terapia_intensiva + ...
            dataset.isolamento_domiciliare)/Npop;
R_uf = (dataset.guariti + dataset.deceduti)/Npop;
S_uf = 1 - I_uf - R_uf;

%% Smoothing and filtering 
% Cubic smoothing 
p = 0.4;
xi = 0:398;

hosp = csaps(xi,log(dataset.ricoverati),p,xi);
ICU = csaps(xi,log(dataset.terapia_intensiva),p,xi);
pos = csaps(xi,log(dataset.isolamento_domiciliare),p,xi);
healed = csaps(xi,log(dataset.guariti),p,xi);
dec = csaps(xi,log(dataset.deceduti),p,xi);

% Apply a Savitzky - Golay filter on the cubic data
order = 5;
framelen = 27;

hosp = sgolayfilt(hosp,order,framelen);
ICU = sgolayfilt(ICU,order,framelen);
pos = sgolayfilt(pos,order,framelen);
healed = sgolayfilt(healed,order,framelen);
dec = sgolayfilt(dec,order,framelen);

hosp = exp(hosp);
ICU = exp(ICU);
pos = exp(pos);
healed = exp(healed);
dec = exp(dec);

t1 = datetime(2020,8,31);
t2 = datetime(2021,10,3);

t = (t1:caldays(1):t2)';

I = (hosp+ICU+pos)'/Npop;
R = (healed+dec)'/Npop;
S = 1 - I - R;
%%

mean_S = mean(S-S_uf);
std_S = std(S-S_uf);

mean_I = mean(I-I_uf);
std_I = std(I-I_uf);

%% Identifiability
opts = optimoptions("lsqnonlin", "OptimalityTolerance",1e-6);
[opt, V] = lsqnonlin(@(par) error_fcn(1:length(S), [S, I], [0 1], par), [0.1,0.1], [],[],[],[],[],[],[],opts);
alpha = opt(1);
gamma = opt(2);

r = 1e-2:1e-2:1;
DeltaV = zeros(size(r));
for R = r
    [opt_pen, V_pen] = lsqnonlin(@(par) pen_fcn(1:length(S), [S, I], [0 1], opt, 1, par), opt, [],[],[],[],[],[],[],opts);
    DeltaV(R == r) = V_pen-V;
end

%% EKF
x0 = [S(1); I(1); 0.1];
ekf = extendedKalmanFilter(@(x) state_transition_NL(x, gamma), ...
                           @(x) x(1:2), ...
                           x0);

ekf.StateCovariance = diag([0.01, 0.01, 100]);
ekf.ProcessNoise = diag([1e-6, 1e-6, 1e-3]);
ekf.MeasurementNoise = diag([std_S, std_I]);

x = [x0'; zeros(length(t)-1,3)];
x_lb = [x0'; zeros(length(t)-1,3)];
x_ub = [x0'; zeros(length(t)-1,3)];
for i = 2:length(t)
    predict(ekf);
    if (mod(i,7) == 0 || i <= 14)
        correct(ekf,[S(i); I(i)]);
    end
    x(i,:) = ekf.State';
    x_ub(i,:) = ekf.State' + 1.96*ones(1,3)*ekf.StateCovariance;
    x_lb(i,:) = ekf.State' - 1.96*ones(1,3)*ekf.StateCovariance;
end

%% Simulation and plotting

fig1 = figure(1);
fig1.Color = 'w';

subplot(311)
plot(t, [S, x(:,1)])
hold on
patch([t; flipud(t)], [x_lb(:,1); flipud(x_ub(:,1))], plot_colors{2}, 'FaceAlpha', 0.5, 'EdgeColor', 'none')
xlim(t([1,end]))
ylabel('$S(t), \hat S(t)$')
legend('$S$', '$\hat S$')

subplot(312)
plot(t, [I, x(:,2)])
hold on
patch([t; flipud(t)], [x_lb(:,2); flipud(x_ub(:,2))], plot_colors{2}, 'FaceAlpha', 0.5, 'EdgeColor', 'none')
xlim(t([1,end]))
ylabel('$I(t), \hat I(t)$')
legend('$I$', '$\hat I$')

subplot(313)
plot(t, x(:,3), 'Color', plot_colors{2});
hold on
patch([t; flipud(t)], [x_lb(:,3); flipud(x_ub(:,3))], plot_colors{2}, 'FaceAlpha', 0.5, 'EdgeColor', 'none')
ylabel('$\hat \alpha(t)$')
xlim(t([1,end]))

%% External functions
function dx = SIR(par, t, x)
    alpha = par(1);
    delta = par(2);
    dx = [     -alpha*x(1)*x(2); ...
          alpha*x(1)*x(2)-delta*x(2)];
end

function e = error_fcn(data_times,data_values,C, par)
    opts = odeset('RelTol',1e-3,'AbsTol',1e-6);
    [t,x] = ode45(@(t,x) SIR(par,t,x), data_times(2:end), data_values(1,:)', opts);
    sigma = zeros(size(x));
    for x_ = x
        sigma(x_ == x) = sqrt(length(t)*(1e-6+1e-3*x_));
    end
    dx = (x - data_values(2:end,:))./sigma;
    e = zeros(length(t),size(C,1));
    for i = 1:length(t)
        e(i,:) = dx(i,:)*C';
    end
end

function e = pen_fcn(data_times,data_values,C, par_star, R, par)
    opts = odeset('RelTol',1e-3,'AbsTol',1e-6);
    [t,x] = ode45(@(t,x) SIR(par,t,x), data_times(2:end), data_values(1,:)', opts);
    sigma = zeros(size(x));
    for x_ = x
        sigma(x_ == x) = sqrt(length(t)*(1e-6+1e-3*x_));
    end
    dx = (x - data_values(2:end,:))./sigma;
    e = zeros(length(t),size(C,1));
    for i = 1:length(t)
        e(i,:) = dx(i,:)*C';
    end
    e = [e(:); sqrt(sum((par-par_star).^2)/R-1)]; 
end

function x = state_transition_NL(x0, gamma)
    [~,x] = ode45(@(t,x) [ -x(3)*x(1)*x(2); x(3)*x(1)*x(2) - gamma*x(2); 0], [0 1], x0);
    x = x(end,:)';
end