function [Phi] = chernick_J2_stm(singular_oe, t, rP, mu, J2)
    %{
    chernick_J2_stm calculates the J2 STM for ROEs defined in Chernick 2.4 
    
    Inputs:
        - singular_oe: keplerian orbit elements [a, e, i, RAAN, w, M] 
        - t: time at which to evaluate the STM (s)
        - Re: radius of the primary attractor (m)
        - mu: Gravitational parameter of the primary attractor (m^3/s^2)
        - J2: J2 coefficient of the primary attractor
        
    Outputs: 
        - Phi: STM evaluated at t
    %}
end

