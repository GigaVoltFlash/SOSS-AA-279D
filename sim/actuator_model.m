function delta_v_real = actuator_model(delta_v_val, apply_noise)
    NOISE_PERCENTAGE = 10;
    if apply_noise
        if any(abs(delta_v_val) > 1e-6) % more than 1e-6 km/s
            noise = (rand()*2 - 1)*NOISE_PERCENTAGE*0.01 * delta_v_val; % noise in [-2%, +2%] of input
            delta_v_real = delta_v_val + noise;
        else
            delta_v_real = delta_v_val;
        end
    else
        delta_v_real = delta_v_val;
    end
end