function [f_syms, f_hand] = getdynamics(eqns,vars,coefs)
    
    nzeros = length(vars)-length(eqns);
    f_syms = [eqns; zeros(nzeros,1)];

    eqnsStr = arrayfun(@(eq) char(eq), f_syms, 'UniformOutput', false);
    
    % replacing all Vars with x values
    replacedEqnsStr = eqnsStr;
    for i = 1:length(eqnsStr)
        for j = 1:length(vars)
            pattern = ['\<', char(vars(j)), '\>'];
            replacement = ['x(', num2str(j), ')'];
            replacedEqnsStr{i} = regexprep(replacedEqnsStr{i}, pattern, replacement);
        end
    end
    
    % replacing all coefs with x values

    for ii = 1:length(eqnsStr)
        for jj = 1:length(coefs)
            pattern = ['\<', char(coefs(jj)), '\>'];
            replacement = ['par(', num2str(jj), ')'];
            replacedEqnsStr{ii} = regexprep(replacedEqnsStr{ii}, pattern, replacement);
        end
    end

    eqnStr = strjoin(replacedEqnsStr, '; ');
    funcStr = ['@(x,par)[', eqnStr, ']'];
    f_hand = str2func(funcStr);

end