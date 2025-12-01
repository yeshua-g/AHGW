function ds = theis_der_dim(params, t, r)
            T = params(1); S = params(2);
            t_d = t .* T ./ (r.^2 .* S);
            ds = 1/2*exp(-1./(4*t_d));
end
