function c_over_csat = solve_equation(x, t, L, D)

      % Calculate the summation term
    summation = 0;
    
    for n = 1:1000  % Adjust the number of terms as needed
        lambda_n = (2 * n - 1) / 2 * pi; 
        term =  sin(lambda_n) / (2*lambda_n + sin(2*lambda_n)) * cos(lambda_n * x / L) * exp(-D * lambda_n^2 / L^2 * t);
        summation = summation + term;
    end

    % Calculate c(x, t) / c_sat(P_eq)
    c_over_csat = 1 - 4*summation;
    if c_over_csat<0 || c_over_csat==0 
        c_over_csat=1e-16;
    end


end