function [par_opt, DeltaV] = identif_test(dyn, C, t, data, x0, par0, varargin)

    import casadi.*
    opts = optimoptions("lsqnonlin", "OptimalityTolerance",1e-6);
    [par_opt, V] = lsqnonlin(@(par) error_fcn(dyn, t, data, x0, C, par, false), par0, [],[],[],[],[],[],[],opts);
    
    if ~isempty(varargin)
        r = varargin{1};
    else
        r = 1e-2:1e-2:1;
    end
    DeltaV = zeros(size(r));
    for R = r
        [~, V_pen] = lsqnonlin(@(par) error_fcn(dyn, t, data, x0, C, par, true, R, par_opt), par_opt, [],[],[],[],[],[],[],opts);
        DeltaV(R == r) = V_pen-V;
    end
end

