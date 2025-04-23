function geom_vals = integration_constants_to_geometry(K, a)
    % geom_vals = [2*a*K23, a*K56, a*K1]
    geom_vals = zeros(3);
    geom_vals(1) = 2*a*sqrt(K(2)^2 + K(3)^2);
    geom_vals(2) = a*sqrt(K(5)^2 + K(6)^2);
    geom_vals(3) = a*K(1);
end