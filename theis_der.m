function ds = theis_der(params, t, Q, r)
    T = params(1); S = params(2);
    u = (r^2 * S) ./ (4 * T * t);
    ds = (Q ./ (4 * pi * T)) .* exp(-u);
end
