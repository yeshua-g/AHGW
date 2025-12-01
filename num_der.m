function dsdlogt = num_der(t, s)
    % t: time vector (lunghezza n)
    % s: drawdown vector
    n = length(t);
    dsdlogt = zeros(size(s));
    
    %provided formula for middle pints
    for i = 2:n-1
        t_im1 = t(i-1); t_i = t(i); t_ip1 = t(i+1);
        s_im1 = s(i-1); s_i = s(i); s_ip1 = s(i+1);
        
        L_ip1_i   = log(t_ip1 / t_i);
        L_i_im1   = log(t_i / t_im1);
        L_ip1_im1 = log(t_ip1 / t_im1);
        
        term1 = (L_i_im1 * s_ip1) / (L_ip1_i * L_ip1_im1);
        term2 = (log((t_ip1 * t_im1) / (t_i^2)) * s_i) / (L_ip1_i * L_i_im1);
        term3 = (L_ip1_i * s_im1) / (L_i_im1 * L_ip1_im1);
        
        dsdlogt(i) = term1 + term2 - term3;
    end
    
    % first and last point: simple differences
    dsdlogt(1) = (s(2) - s(1)) / (log(t(2)) - log(t(1)));
    dsdlogt(n) = (s(n) - s(n-1)) / (log(t(n)) - log(t(n-1)));
end
