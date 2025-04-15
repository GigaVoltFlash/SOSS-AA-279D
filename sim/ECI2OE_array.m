function [a,e,i,RAAN,omega,nu] = ECI2OE_array(r_ECI,v_ECI)
    N = length(r_ECI);
    a = zeros(1,N);
    e = zeros(1,N);
    i = zeros(1,N);
    RAAN = zeros(1,N); 
    omega = zeros(1,N);
    nu = zeros(1,N);
    
    for k = 1:N  
        r = r_ECI(k,:);
        v = v_ECI(k,:);
        
        [a(k),e(k),i(k),RAAN(k),omega(k),nu(k)] = ECI2OE(r,v);
    end
end