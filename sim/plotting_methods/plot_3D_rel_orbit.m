function [] = plot_3D_rel_orbit(SV2_rel_pos1, SV3_rel_pos1, fig_path)
    x1_1 = SV2_rel_pos1(:,1); y1_1 = SV2_rel_pos1(:,2); z1_1 = SV2_rel_pos1(:,3);
    x2_1 = SV3_rel_pos1(:,1); y2_1 = SV3_rel_pos1(:,2); z2_1 = SV3_rel_pos1(:,3);
    
    figure();
    plot3(x1_1, y1_1, z1_1, 'b');
    hold on;


    xlabel('Radial [km]');
    ylabel('Tangential [km]');
    zlabel('Normal [km]');
    title('SV2 3D Relative Orbits with HCW')
    grid on;
    axis equal;
    %legend('YA Dif. Eq.', 'YA Mapping');
    hold off;

    figure();
    plot3(x2_1, y2_1, z2_1, 'b');
    hold on;


    xlabel('Radial [km]');
    ylabel('Tangential [km]');
    zlabel('Normal [km]');
    title('SV3 3D Relative Orbits with HCW')
    grid on;
    axis equal;
    %legend('YA Dif. Eq.', 'YA Mapping');
    hold off;

    saveas(gcf, fig_path);
end