function [a_result, e_x, e_y, i_result, RAAN_result, u] = OE2quasi_nonsing_array(a,e,i,RAAN,w,M)
    % Takes in deg, outputs in deg

    N = length(a);
    a_result = zeros(1,N);
    e_x = zeros(1,N);
    e_y = zeros(1,N);
    i_result = zeros(1,N);
    RAAN_result = zeros(1,N);
    u = zeros(1,N);

    for k = 1:N
        [a_result(k),e_x(k),e_y(k),i_result(k),RAAN_result(k),u(k)] = OE2quasi_nonsing(a(k),e(k),i(k),RAAN(k),w(k),M(k));
    end

end