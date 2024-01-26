clc; clear; close all;
default_paper;

%% Dynamics
a = sym('a', 'real');
c = sym('c', 'real');
S = sym('S', 'real');
I = sym('I', 'real');
R = sym('R', 'real');
var = {S; I; a; c};
n = length(var);

dyn = { -a*S*I      ; ...
        a*S*I - c*I ; ...
        0           ; ...
        0           };

C = [1, 1, 0, 0];

%% Observability
J = sym('x', n);
O = sym('x', n);

for i = 1:n
    for j = 1:n
        J(i,j) = diff(dyn{i},var{j});
    end
end

O(1,:) = C;
for i = 2:n
    O(i,:) = O(i-1,:)*J;
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
T = [o,uo];
simplify(C*T)
simplify(T*[S;I;a;c])
simplify(T\J*T)