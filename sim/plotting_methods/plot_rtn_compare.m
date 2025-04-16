% Assumes:
% - r_RTN: [nx3] position in RTN frame
% - v_RTN: [nx3] velocity in RTN frame
% - t:     [nx1] time vector

function [] = plot_rtn_compare(t, t_orbit, r_RTN_1, r_RTN_2, v_RTN_1, v_RTN_2, label1, label2, title_string, fig_path)
    figure('Units', 'normalized', 'OuterPosition', [0 0 1 1]);  % Fullscreen
    
    titles = {'Radial (R)', 'Transverse (T)', 'Normal (N)'};
    
    for i = 1:3
        % Position subplot (left column)
        subplot(3, 2, 2*i - 1);
        plot(t/t_orbit, r_RTN_1(:, i) - r_RTN_2(:, i), 'b', 'LineWidth', 1.5);
        ylabel(['r_{RTN} - ' titles{i} ' [km]']);
        if i == 1
            title('Error in RTN Position Components');
        end
        if i == 3
            xlabel('Time [s]');
        end
        grid on;
    
        % Velocity subplot (right column)
        subplot(3, 2, 2*i);
        plot(t/t_orbit, v_RTN_1(:, i) - v_RTN_2(:, i), 'r', 'LineWidth', 1.5);
        ylabel(['v_{RTN} - ' titles{i} ' [km/s]']);
        if i == 1
            title('Error in RTN Velocity Components');
        end
        if i == 3
            xlabel('Time [s]');
        end
        grid on;
    end
    
    % xlabel('Time [s]');
    sgtitle(title_string);
    saveas(gcf, fig_path);

    %%% Additional code if you want to plot both the actual values instead
    %%% of just the error.
    % figure('Units', 'normalized', 'OuterPosition', [0 0 1 1]);  % Fullscreen
    % 
    % titles = {'Radial (R)', 'Transverse (T)', 'Normal (N)'};
    % 
    % for i = 1:3
    %     % Position subplot (left column)
    %     subplot(3, 2, 2*i - 1);
    %     plot(t/t_orbit, r_RTN_1(:, i), 'b', 'LineWidth', 1.5); hold on;
    %     plot(t/t_orbit, r_RTN_2(:, i), 'r--', 'LineWidth', 1.5);
    %     ylabel(['r_{RTN} - ' titles{i} ' [km]']);
    %     if i == 1
    %         title('RTN Position Components');
    %         legend(label1, label2);
    %     end
    %     grid on;
    % 
    %     % Velocity subplot (right column)
    %     subplot(3, 2, 2*i);
    %     plot(t/t_orbit, v_RTN_1(:, i), 'b', 'LineWidth', 1.5); hold on;
    %     plot(t/t_orbit, v_RTN_2(:, i), 'r--', 'LineWidth', 1.5);
    %     ylabel(['v_{RTN} - ' titles{i} ' [km/s]']);
    %     if i == 1
    %         title('RTN Velocity Components');
    %         legend(label1, label2);
    %     end
    %     grid on;
    % end
    % 
    % xlabel('Time [s]');
    % sgtitle(title_string);
