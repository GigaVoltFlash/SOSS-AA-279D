function [a,e,i,RAAN,omega,nu,M] = ECI2OE_array(r_ECI,v_ECI)
    if size(r_ECI, 1) == 1  % Single vector input
        [a,e,i,RAAN,omega,nu,M] = ECI2OE(r_ECI, v_ECI);
    else                    % Multiple vectors
        N = size(r_ECI, 1);
        a = zeros(1,N);
        e = zeros(1,N);
        i = zeros(1,N);
        RAAN = zeros(1,N); 
        omega = zeros(1,N);
        nu = zeros(1,N);
        M = zeros(1,N); 

        for k = 1:N  
            r = r_ECI(k,:);
            v = v_ECI(k,:);
            [a(k),e(k),i(k),RAAN(k),omega(k),nu(k),M(k)] = ECI2OE(r,v);
        end
    end
end