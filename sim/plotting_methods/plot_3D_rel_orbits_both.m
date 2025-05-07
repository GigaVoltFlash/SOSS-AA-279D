function [] = plot_3D_rel_orbits_both(SV2_rel_pos1, SV2_rel_pos2, SV3_rel_pos1, SV3_rel_pos2,fig_path)
    x1_1 = SV2_rel_pos1(:,1); y1_1 = SV2_rel_pos1(:,2); z1_1 = SV2_rel_pos1(:,3);
    x2_1 = SV3_rel_pos1(:,1); y2_1 = SV3_rel_pos1(:,2); z2_1 = SV3_rel_pos1(:,3);
    x1_2 = SV2_rel_pos2(:,1); y1_2 = SV2_rel_pos2(:,2); z1_2 = SV2_rel_pos2(:,3);
    x2_2 = SV3_rel_pos2(:,1); y2_2 = SV3_rel_pos2(:,2); z2_2 = SV3_rel_pos2(:,3);
    
    figure();
    plot3(x1_1, y1_1, z1_1, 'b');
    hold on;
    plot3(x1_2, y1_2, z1_2, 'r');

    xlabel('Radial [km]');
    ylabel('Tangential [km]');
    zlabel('Normal [km]');
    title('SV2 3D Relative Orbits')
    grid on;
    axis equal;
    legend('YA Dif. Eq.', 'YA Mapping');
    hold off;

    figure();
    plot3(x2_1, y2_1, z2_1, 'b');
    hold on;
    plot3(x2_2, y2_2, z2_2, 'r');

    xlabel('Radial [km]');
    ylabel('Tangential [km]');
    zlabel('Normal [km]');
    title('SV3 3D Relative Orbits')
    grid on;
    axis equal;
    legend('YA Dif. Eq.', 'YA Mapping');
    hold off;

    saveas(gcf, fig_path);
end