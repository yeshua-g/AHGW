function s = theis_dimensionsless(params,t,r,r_w)
% Se r_w non è stato passato, lo imposto uguale a r
    if nargin < 4
        r_w = r;
    end
    T = params(1); S = params(2);
    r_d = r./r_w;
    t_d = t .* T ./ (r_w.^2 .* S);
    s = 1/2*expint(r_d.^2./(4*t_d));
end

