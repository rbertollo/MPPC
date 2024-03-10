function mape = mapeFunc(y_true, y_pred)
    % check on input vector size
    if numel(y_true) ~= numel(y_pred)
        error('Input vectors must have the same size.');
    end
    
    
    abs_perc_err = abs((y_true - y_pred) ./ y_true);
    abs_perc_err(isnan(abs_perc_err)) = 0;
    
    % Calculate MAPE
    mape = mean(abs_perc_err) * 100;
end
