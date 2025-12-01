function s = theis_model_rec(params, Q, tr, tp)
    T = params;
    s = 2.3*(Q ./ (4 * pi * T)) .* log10((tr+tp)./tr);
end

