function x = state_transition(f_ODE,x0)
    [~,x] = ode45(f_ODE, [0 1], x0);
    x = x(end,:)';
end