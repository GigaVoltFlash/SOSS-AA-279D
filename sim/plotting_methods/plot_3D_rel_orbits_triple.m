function [] = plot_3D_rel_orbits_triple(SV2_rel_pos1, SV2_rel_pos2, SV2_rel_pos3, SV3_rel_pos1, SV3_rel_pos2, SV3_rel_pos3,title_string,fig_path)
    x1_1 = SV2_rel_pos1(:,1); y1_1 = SV2_rel_pos1(:,2); z1_1 = SV2_rel_pos1(:,3);
    x2_1 = SV3_rel_pos1(:,1); y2_1 = SV3_rel_pos1(:,2); z2_1 = SV3_rel_pos1(:,3);
    x1_2 = SV2_rel_pos2(:,1); y1_2 = SV2_rel_pos2(:,2); z1_2 = SV2_rel_pos2(:,3);
    x2_2 = SV3_rel_pos2(:,1); y2_2 = SV3_rel_pos2(:,2); z2_2 = SV3_rel_pos2(:,3);
    x1_3 = SV2_rel_pos3(:,1); y1_3 = SV2_rel_pos3(:,2); z1_3 = SV2_rel_pos3(:,3);
    x2_3 = SV3_rel_pos3(:,1); y2_3 = SV3_rel_pos3(:,2); z2_3 = SV3_rel_pos3(:,3);
    
    figure();

    t = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(t, title_string, 'FontWeight', 'bold');

    ax1 = nexttile;
    plot3(x1_1, y1_1, z1_1, 'b');
    hold on;
    plot3(x1_2, y1_2, z1_2, 'r');
    plot3(x1_3, y1_3, z1_3, 'g');
    xlabel('Radial [km]');
    ylabel('Tangential [km]');
    zlabel('Normal [km]');
    title('SV2 3D Relative Orbits')
    grid on;
    axis equal;
    %legend('YA Dif. Eq.', 'YA Mapping', 'YA FERM');
    %legend({'YA Dif. Eq.', 'YA Mapping', 'YA FERM'}, 'Location', 'northeastoutside');
    hold off;

    ax2 = nexttile;
    plot3(x2_1, y2_1, z2_1, 'b');
    hold on;
    plot3(x2_2, y2_2, z2_2, 'r');
    plot3(x2_3, y2_3, z2_3, 'g');
    xlabel('Radial [km]');
    ylabel('Tangential [km]');
    zlabel('Normal [km]');
    title('SV3 3D Relative Orbits')
    grid on;
    axis equal;
    %legend('YA Dif. Eq.', 'YA Mapping', 'YA FERM');
    %legend({'YA Dif. Eq.', 'YA Mapping', 'YA FERM'}, 'Location', 'northeastoutside');
    hold off;

    lgd = legend(ax2, {'YA Dif. Eq.', 'YA Mapping', 'YA FERM'}, 'Location', 'southoutside', 'Orientation', 'horizontal');
    lgd.Layout.Tile = 'south';  % Place legend under both plots
    % Save figure
    saveas(gcf, fig_path);
end