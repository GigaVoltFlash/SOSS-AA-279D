function [a_2,e_2,i_2,RAAN_2,omega_2,nu_2] = compute_and_plot_orbit_params(r_ECI_no_j2,v_ECI_no_j2,r_ECI_with_j2,v_ECI_with_j2,t_1,t_2)
    N = length(r_ECI_no_j2);
    a_1 = zeros(1,N);
    e_1 = zeros(1,N);
    i_1 = zeros(1,N);
    RAAN_1 = zeros(1,N); 
    omega_1 = zeros(1,N);
    nu_1 = zeros(1,N);
    h_1 = zeros(3,N);
    e_vec_1 = zeros(3,N);
    epsilon_1 = zeros(1,N);
    a_2 = zeros(1,N);
    e_2 = zeros(1,N);
    i_2 = zeros(1,N);
    RAAN_2 = zeros(1,N); 
    omega_2 = zeros(1,N);
    nu_2 = zeros(1,N);
    h_2 = zeros(3,N);
    e_vec_2 = zeros(3,N);
    epsilon_2 = zeros(1,N);
    
    for k = 1:N
        r_1 = r_ECI_no_j2(k,:);
        v_1 = v_ECI_no_j2(k,:);
        r_2 = r_ECI_with_j2(k,:);
        v_2 = v_ECI_with_j2(k,:);
    
        [a_1(k),e_1(k),i_1(k),RAAN_1(k),omega_1(k),nu_1(k),h_1(:,k),e_vec_1(:,k),epsilon_1(k)] = orbit_params(r_1,v_1);
        [a_2(k),e_2(k),i_2(k),RAAN_2(k),omega_2(k),nu_2(k),h_2(:,k),e_vec_2(:,k),epsilon_2(k)] = orbit_params(r_2,v_2);
        
        
    end

    figure;

    subplot(3,2,1)
    plot(t_1, a_1, 'b', t_2, a_2, 'r--')
    xlabel('Time [s]'); ylabel('Semi-major axis [km]')
    title('Semi-major axis'); legend('Without J2','With J2'); grid on
    
    subplot(3,2,2)
    plot(t_1, e_1, 'b', t_2, e_2, 'r--')
    xlabel('Time [s]'); ylabel('Eccentricity')
    title('Eccentricity'); legend('Without J2','With J2'); grid on
    
    subplot(3,2,3)
    plot(t_1, i_1, 'b', t_2, i_2, 'r--')
    xlabel('Time [s]'); ylabel('Inclination [deg]')
    %ylim([95 105])
    title('Inclination'); legend('Without J2','With J2'); grid on
    
    subplot(3,2,4)
    plot(t_1, RAAN_1, 'b', t_2, RAAN_2, 'r--')
    xlabel('Time [s]'); ylabel('RAAN [deg]')
    title('RAAN'); legend('Without J2','With J2'); grid on
    
    subplot(3,2,5)
    plot(t_1, omega_1, 'b', t_2, omega_2, 'r--')
    xlabel('Time [s]'); ylabel('Argument of Perigee [deg]')
    title('Argument of Perigee'); legend('Without J2','With J2'); grid on
    
    subplot(3,2,6)
    plot(t_1, nu_1, 'b', t_2, nu_2, 'r--')
    xlabel('Time [s]'); ylabel('True Anomaly [deg]')
    title('True Anomaly'); legend('Without J2','With J2'); grid on
    
    %sgtitle('Orbital Elements Comparison (Without J2 vs With J2)')
    
    figure()
    plot(t_1, epsilon_1, 'b', t_2, epsilon_2, 'r--')
    xlabel('Time [s]'); ylabel('\epsilon [km^2/s^2]')
    ylim([-30 -28])
    title('Specific Mechanical Energy'); legend('Without J2','With J2'); grid on
    
    figure;
    subplot(3,1,1)
    plot(t_1, h_1(1,:), 'b', t_2, h_2(1,:), 'r--')
    xlabel('Time [s]'); ylabel('h_x [km^2/s]')
    title('Angular Momentum: x-component'); legend('Without J2','With J2'); grid on
    
    subplot(3,1,2)
    plot(t_1, h_1(2,:), 'b', t_2, h_2(2,:), 'r--')
    xlabel('Time [s]'); ylabel('h_y [km^2/s]')
    title('Angular Momentum: y-component'); legend('Without J2','With J2'); grid on
    
    subplot(3,1,3)
    plot(t_1, h_1(3,:), 'b', t_2, h_2(3,:), 'r--')
    xlabel('Time [s]'); ylabel('h_z [km^2/s]')
    ylim([-9000 -8000])
    title('Angular Momentum: z-component'); legend('Without J2','With J2'); grid on
    
    %sgtitle('Angular Momentum Vector Components Over Time')
    
    figure;
    subplot(3,1,1)
    plot(t_1, e_vec_1(1,:), 'b', t_2, e_vec_2(1,:), 'r--')
    xlabel('Time [s]'); ylabel('e_x')
    title('Eccentricity Vector: x-component'); legend('Without J2','With J2'); grid on
    
    subplot(3,1,2)
    plot(t_1, e_vec_1(2,:), 'b', t_2, e_vec_2(2,:), 'r--')
    xlabel('Time [s]'); ylabel('e_y')
    title('Eccentricity Vector: y-component'); legend('Without J2','With J2'); grid on
    
    subplot(3,1,3)
    plot(t_1, e_vec_1(3,:), 'b', t_2, e_vec_2(3,:), 'r--')
    xlabel('Time [s]'); ylabel('e_z')
    title('Eccentricity Vector: z-component'); legend('Without J2','With J2'); grid on
    
    %sgtitle('Eccentricity Vector Components Over Time')
end
