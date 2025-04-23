function [] = plot_3D_rel_orbit(SV2_rel_pos1, SV3_rel_pos1, fig_path_SV2, fig_path_SV3)
    x1 = SV2_rel_pos1(:,1); y1 = SV2_rel_pos1(:,2); z1 = SV2_rel_pos1(:,3);
    x2 = SV3_rel_pos1(:,1); y2 = SV3_rel_pos1(:,2); z2 = SV3_rel_pos1(:,3);
    
    figure();
    plot3(x1, y1, z1, 'b');
    hold on;


    xlabel('Radial [km]');
    ylabel('Tangential [km]');
    zlabel('Normal [km]');
    title('SV2 3D Relative Orbits with HCW')
    grid on;
    axis equal;
    %legend('YA Dif. Eq.', 'YA Mapping');
    hold off;
    saveas(gcf, fig_path_SV2);

    figure();
    plot3(x2, y2, z2, 'b');
    hold on;


    xlabel('Radial [km]');
    ylabel('Tangential [km]');
    zlabel('Normal [km]');
    title('SV3 3D Relative Orbits with HCW')
    grid on;
    axis equal;
    %legend('YA Dif. Eq.', 'YA Mapping');
    hold off;

    saveas(gcf, fig_path_SV3);
end