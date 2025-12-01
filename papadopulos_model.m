

function s = papadopulos_model(p, t, Q, r)
    T = p(1); S = p(2); C = p(3);
    % Dimensionless parameters
    u = (r^2 * S) ./ (4 * T * t);
    % Theis well function
    W = expint(u);
    % Papadopulos-Cooper correction term
    % Dimensionless storage parameter
    sigma = (C * T) / (r^2 * S);
    % Compute correction using series approximation
    % s = (Q/(4*pi*T)) * (W + correction)
    correction = zeros(size(t));
    for i = 1:length(t)
        % Integral approximation for correction
        % Papadopulos-Cooper uses Laplace inversion; here we use a series
        % sum_{n=1}^{N} exp(-n^2 * pi^2 * t(i) / sigma)
        N = 50; % number of terms for accuracy
        sum_term = 0;
        for n = 1:N
            sum_term = sum_term + exp(-(n^2 * pi^2 * t(i)) / sigma);
        end
        correction(i) = (2 / pi) * sum_term;
    end
    s = (Q / (4 * pi * T)) .* (W + correction);
end
