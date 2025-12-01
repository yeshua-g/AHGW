function s = theis_model(params, t, Q, r)
    T = params(1); S = params(2);
    u = (r^2 * S) ./ (4 * T * t);
    s = (Q ./ (4 * pi * T)) .* expint(u);
end

