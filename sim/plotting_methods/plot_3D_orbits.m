function [] = plot_3D_orbits(r_ECI_1, r_ECI_2, title_1, title_2)
    x1 = r_ECI_1(:,1);
    y1 = r_ECI_1(:,2);
    z1 = r_ECI_1(:,3);
    x2 = r_ECI_2(:,1);
    y2 = r_ECI_2(:,2);
    z2 = r_ECI_2(:,3);
    
    figure();
    plot3(x1, y1, z1, 'b');
    hold on;
    plot3(x2, y2, z2, 'r');
    scatter3(0, 0, 0, 100, 'go', 'MarkerFaceColor', 'r');  % Central body
    scatter3(x1(1), y1(1), z1(1), 70, 'bo', 'MarkerFaceColor', 'b');  % Start
    scatter3(x1(end), y1(end), z1(end), 70, 'mo', 'MarkerFaceColor', 'm');  % End
    scatter3(x2(1), y2(1), z2(1), 70, 'ro', 'MarkerFaceColor', 'r');  % Start
    scatter3(x2(end), y2(end), z2(end), 70, 'yo', 'MarkerFaceColor', 'y');  % End
    xlabel('x');
    ylabel('y');
    zlabel('z');
    title('Chief and Deputy Orbits')
    grid on;
    axis equal;
    legend(title_1, title_2,'Central Body', 'Start Point', 'End Point');
    hold off;
end