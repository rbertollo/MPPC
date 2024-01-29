%% Observability function in SYSM
% Fucntion that calculates the Observabillity in symbolic given:
% J = Symbolic Jacobian of the system
% C = Meaurement matrix

function O = obsv_sym(J,C) 
    n = size(J,1);
    O = C;
    for i = 1:n-1
        newrow = C * J^i;
        O = [O; newrow];
    end
end