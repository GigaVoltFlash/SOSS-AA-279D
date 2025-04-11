function [] = plot_osc_mean_oe(a_2,e_2,i_2,RAAN_2,omega_2,nu_2,t_2,a_3,e_3,i_3,RAAN_3,omega_3,nu_3,t_3)
    figure;

    subplot(3,2,1)
    plot(t_2, a_2, 'b', t_3, a_3, 'r--')
    xlabel('Time [s]'); ylabel('Semi-major axis [km]')
    title('Semi-major axis'); legend('Osculating with J2','Mean with J2'); grid on
    
    subplot(3,2,2)
    plot(t_2, e_2, 'b', t_3, e_3, 'r--')
    xlabel('Time [s]'); ylabel('Eccentricity')
    title('Eccentricity'); legend('Osculating with J2','Mean with J2'); grid on
    
    subplot(3,2,3)
    plot(t_2, i_2, 'b', t_3, i_3, 'r--')
    xlabel('Time [s]'); ylabel('Inclination [deg]')
    title('Inclination'); legend('Osculating with J2','Mean with J2'); grid on
    
    subplot(3,2,4)
    plot(t_2, RAAN_2, 'b', t_3, RAAN_3, 'r--')
    xlabel('Time [s]'); ylabel('RAAN [deg]')
    title('RAAN'); legend('Osculating with J2','Mean with J2'); grid on
    
    subplot(3,2,5)
    plot(t_2, omega_2, 'b', t_3, omega_3, 'r--')
    xlabel('Time [s]'); ylabel('Argument of Perigee [deg]')
    title('Argument of Perigee'); legend('Osculating with J2','Mean with J2'); grid on
    
    subplot(3,2,6)
    plot(t_2, nu_2, 'b', t_3, nu_3, 'r--')
    xlabel('Time [s]'); ylabel('True Anomaly [deg]')
    title('True Anomaly'); legend('Osculating with J2','Mean with J2'); grid on
    
    %sgtitle('Orbital Elements Comparison with J2 (Osculating vs Mean)')

end