function newFunction = ParamEvalFunc(funcHandle, values)
    
    funcStr = func2str(funcHandle);
    
    for i = 1:length(values)
        pattern = ['par(', num2str(i), ')'];
        replacement = num2str(values(i));
        funcStr = strrep(funcStr, pattern, replacement);
    end
    newFunction = str2func(funcStr);
end