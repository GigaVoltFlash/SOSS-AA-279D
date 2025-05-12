function create_roe_3d_animation(SV2_rel_pos, SV3_rel_pos, title_string)
    x1_1 = SV2_rel_pos(:,1); y1_1 = SV2_rel_pos(:,2); z1_1 = SV2_rel_pos(:,3);
    x2_1 = SV3_rel_pos(:,1); y2_1 = SV3_rel_pos(:,2); z2_1 = SV3_rel_pos(:,3);
    
     % === Setup ===
    fig = figure();
    hold on; grid on; axis equal;
    
    xlabel('Radial [km]');
    ylabel('Tangential [km]');
    zlabel('Normal [km]');
    sgtitle(title_string, 'FontWeight', 'bold');

    view(3);
    
    % Chief (SV1) at origin
    plot3(0, 0, 0, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6, 'DisplayName', 'SV1 / Chief');
    
    % Initialize trajectory lines and moving markers
    h_sv2 = plot3(NaN, NaN, NaN, 'b-', 'LineWidth', 1.5, 'DisplayName', 'SV2 / Watcher');
    h_sv3 = plot3(NaN, NaN, NaN, 'r-', 'LineWidth', 1.5, 'DisplayName', 'SV3 / Docker');
    
    pt_sv2 = plot3(NaN, NaN, NaN, 'bo', 'MarkerSize', 6, 'MarkerFaceColor', 'b');
    pt_sv3 = plot3(NaN, NaN, NaN, 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r');
    
    legend('Location', 'best');
    
    % === Video Recording Setup ===
    video_filename = 'figures/PS5/relative_orbit_classic.mp4';
    v = VideoWriter(video_filename, 'MPEG-4');
    v.FrameRate = 30;
    open(v);
    
    % === Animation Loop ===
    step = 20;  % Skip 5 time steps per frame
    N = length(x1_1);
    for k = 1:step:N
        % Update trajectories
        set(h_sv2, 'XData', x1_1(1:k), 'YData', y1_1(1:k), 'ZData', z1_1(1:k));
        set(h_sv3, 'XData', x2_1(1:k), 'YData', y2_1(1:k), 'ZData', z2_1(1:k));
    
        % Update moving markers
        set(pt_sv2, 'XData', x1_1(k), 'YData', y1_1(k), 'ZData', z1_1(k));
        set(pt_sv3, 'XData', x2_1(k), 'YData', y2_1(k), 'ZData', z2_1(k));
    
        drawnow;
        frame = getframe(fig);
        writeVideo(v, frame);
    end
    
    close(v);
    disp(['MP4 saved to: ', video_filename]);
end