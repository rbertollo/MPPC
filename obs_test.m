clc; 
clear all ; 
close all ;

%% Dynamics

% a = sym('a', 'real');
% c = sym('c', 'real');
% S = sym('S', 'real');
% I = sym('I', 'real');
% R = sym('R', 'real');
% var = {S; I; a; c};
% n = length(var);
% 
% dyn = { -a*S*I      ; ...
%         a*S*I - c*I ; ...
%         0           ; ...
%         0           };
% 
% C = [1, 1, 0, 0];

%% Dynamics

alpha = sym('alpha', 'real');
gamma = sym('gamma', 'real');
delta = sym('delta', 'real');
sigma = sym('sigma', 'real');
tau = sym('tau', 'real');
lambda = sym('lambda', 'real');
S = sym('S', 'real');
I = sym('I', 'real');
D = sym('D', 'real');
T = sym('T', 'real');
H = sym('H', 'real');
E = sym('E', 'real');

var = {S; I; D; T; H; alpha; gamma; delta};
n = length(var);

dyn = { -S * (alpha * I); ...
         S * (alpha * I) - gamma * I; ...
         I * gamma - D * (lambda + delta); ...
         D * delta - (sigma + tau) * T; ...
        (I + D) * lambda + T * sigma; ...
        % E is removed 
        0           ; ...
        0           ; ...
        0 };

C = [eye(5), zeros(5,3)];
n_meas = size(C,1);

%% Observability
J = sym('x', n);
O = sym('x', n);

for i = 1:n
    for j = 1:n
        J(i,j) = diff(dyn{i},var{j});
    end
end

O(1:n_meas,:) = C;
for i = 1:(n-1)
    O((i*n_meas+1):(i+1)*n_meas,:) = O(((i-1)*n_meas+1):i*n_meas,:)*J;
end
rank(O)

%% Unobservable subspace
syms x [n,1];
syms uo [n,1];
sol = solve([O*x == 0; x1 == 1], x, 'Real', true);
for i = 1:n
    uo(i) = sol.(strcat('x', num2str(i)));
end

%% Observable subspace
syms x [n, 1];
syms o [n, rank(O)];
for i = 1:rank(O)
    sol_found = false;
    j = 1;
    while ~sol_found
        sol = solve([[uo,o(:,1:i-1)]'*x == 0; x(j) == 1], x, 'Real', true);
        if ~isempty(sol.(strcat('x', num2str(i))))
            sol_found = true;
            for k = 1:n
                o(k,i) = sol.(strcat('x', num2str(k)));
            end
            break
        else
            j = j+1;
        end
    end
end

%% "Kalman" decomposition
% T = [o,uo];
% simplify(C*T)
% simplify(T*[var{:}].')
% simplify(T\J*T)