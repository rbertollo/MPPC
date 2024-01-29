function e = error_fcn(dyn, data_times, data_values, x0, C, par)
    
    % Simulate the system trajectories
    opts = odeset('RelTol',1e-3,'AbsTol',1e-6);
    [t,x] = ode45(@(t,x) dyn(x,par), data_times(2:end), x0, opts);
    
    % Compute the measured output
    y = zeros(length(t),size(C,1));
    for i = 1:length(t)
        y(i,:) = x(i,:)*C';
    end
    
    % Compute the std of each point (numerical tolerance of ODE integrator)
    sigma = zeros(size(y));
    for y_ = y
        sigma(y_ == y) = sqrt(length(t)*(1e-6+1e-3*y_));
    end
    
    % Evaluate the difference between simulated points and data points
    dy = (y - data_values(2:end,:))./sigma;
    e = dy(:);
end