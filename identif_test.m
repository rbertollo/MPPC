function [opt, DeltaV] = identif_test(dyn, C, t, data, x0, par0, varargin)
    opts = optimoptions("lsqnonlin", "OptimalityTolerance",1e-6);
    [opt, V] = lsqnonlin(@(par) error_fcn(dyn, t, data, x0, C, par), par0, [],[],[],[],[],[],[],opts);
    
    if nargin == 5
        r = varargin{1};
    else
        r = 1e-2:1e-2:1;
    end
    DeltaV = zeros(size(r));
    for R = r
        [~, V_pen] = lsqnonlin(@(par) pen_fcn(dyn, t, data, x0, C, opt, R, par), opt, [],[],[],[],[],[],[],opts);
        DeltaV(R == r) = V_pen-V;
    end
end

