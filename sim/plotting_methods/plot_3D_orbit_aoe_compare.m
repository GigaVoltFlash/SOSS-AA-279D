function [] = plot_3D_orbit_aoe_compare(r_ECI_1, r_ECI_2, title_1, title_2)
    x1 = r_ECI_1(:,1);
    y1 = r_ECI_1(:,2);
    z1 = r_ECI_1(:,3);
    x2 = r_ECI_2(:,1);
    y2 = r_ECI_2(:,2);
    z2 = r_ECI_2(:,3);
    
    figure('Units', 'normalized', 'OuterPosition', [0 0 1 1]);

    % --- First subplot ---
    ax1 = subplot(1, 2, 1);  % 1 row, 2 columns, 1st plot
    plot3(x1, y1, z1, 'b');
    hold on;
    scatter3(0, 0, 0, 100, 'go', 'MarkerFaceColor', 'r');  % Central body
    scatter3(x1(1), y1(1), z1(1), 70, 'ro', 'MarkerFaceColor', 'g');  % Start
    scatter3(x1(end), y1(end), z1(end), 70, 'mo', 'MarkerFaceColor', 'm');  % End
    xlabel('x');
    ylabel('y');
    zlabel('z');
    title(title_1);
    grid on;
    axis equal;
    legend('Orbit', 'Central Body', 'Start Point', 'End Point');
    hold off;

    % --- Second subplot ---
    ax2 = subplot(1, 2, 2);  % 1 row, 2 columns, 2nd plot
    plot3(x2, y2, z2, 'r'); 
    hold on;
    scatter3(0, 0, 0, 100, 'go', 'MarkerFaceColor', 'r');  % Central body
    scatter3(x2(1), y2(1), z2(1), 70, 'ro', 'MarkerFaceColor', 'g');  % Start
    scatter3(x2(end), y2(end), z2(end), 70, 'mo', 'MarkerFaceColor', 'm');  % End
    xlabel('x');
    ylabel('y');
    zlabel('z');
    title(title_2);
    grid on;
    axis equal;
    legend('Orbit', 'Central Body', 'Start Point', 'End Point');
    hold off;
end