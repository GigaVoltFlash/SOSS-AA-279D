% function [] = plot_rtn_compare(t, r_RTN_1, r_RTN_2, v_RTN_1, v_RTN_2, label1, label2, title)
%     figure('Units', 'normalized', 'OuterPosition', [0 0 1 1]);
% 
%     titles = {'Radial (R)', 'Transverse (T)', 'Normal (N)'};
% 
%     % Position components
%     subplot(3, 1, 1);
%     plot(t, v_RTN_1(:, 1), 'b', 'LineWidth', 1.5); hold on;
%     plot(t, v_RTN_2(:, 1), 'r--', 'LineWidth', 1.5);
%     ylabel(['r_{RTN} - ' titles{1} ' [km]']);
%     legend(label1, label2);
%     grid on;
% 
%     % Velocity components
%     subplot(3, 1, 2);
%     plot(t, v_RTN_1(:, 2), 'b', 'LineWidth', 1.5); hold on;
%     plot(t, v_RTN_2(:, 2), 'r--', 'LineWidth', 1.5);
%     ylabel(['v_{RTN} - ' titles{1} ' [km/s]']);
%     grid on;
% 
%     subplot(3, 1, 3);
%     plot(t, v_RTN_1(:, 3), 'b', 'LineWidth', 1.5); hold on;
%     plot(t, v_RTN_2(:, 3), 'r--', 'LineWidth', 1.5);
%     ylabel(['v_{RTN} - ' titles{2} ' [km/s]']);
%     grid on;
% 
%     xlabel('Time [s]');
%     sgtitle(title);
% 
%     figure('Units', 'normalized', 'OuterPosition', [0 0 1 1]);
% 
% 
%     % Position components
%     subplot(3, 1, 1);
%     plot(t, r_RTN_1(:, 1) - r_RTN_2(:, 1), 'b', 'LineWidth', 1.5); hold on;
%     ylabel(['Error in r_{RTN} - ' titles{1} ' [km]']);
%     grid on;
% 
%     % Velocity components
%     subplot(3, 1, 2);
%     plot(t, v_RTN_1(:, 1) - v_RTN_2(:, 1), 'b', 'LineWidth', 1.5); hold on;
%     ylabel(['Error in v_{RTN} - ' titles{1} ' [km/s]']);
%     grid on;
% 
%     subplot(3, 1, 3);
%     plot(t, v_RTN_1(:, 2) - v_RTN_2(:, 2), 'b', 'LineWidth', 1.5); hold on;
%     ylabel(['Error in v_{RTN} - ' titles{2} ' [km/s]']);
%     grid on;
% 
%     xlabel('Time [s]');
%     sgtitle(title);
% end

% Assumes:
% - r_RTN: [nx3] position in RTN frame
% - v_RTN: [nx3] velocity in RTN frame
% - t:     [nx1] time vector

function [] = plot_rtn_compare(t, r_RTN_1, r_RTN_2, v_RTN_1, v_RTN_2, label1, label2, title_string)
    figure('Units', 'normalized', 'OuterPosition', [0 0 1 1]);  % Fullscreen
    
    titles = {'Radial (R)', 'Transverse (T)', 'Normal (N)'};
    
    for i = 1:3
        % Position subplot (left column)
        subplot(3, 2, 2*i - 1);
        plot(t, r_RTN_1(:, i) - r_RTN_2(:, i), 'b', 'LineWidth', 1.5);
        ylabel(['r_{RTN} - ' titles{i} ' [km]']);
        if i == 1
            title('Error in RTN Position Components');
        end
        grid on;
    
        % Velocity subplot (right column)
        subplot(3, 2, 2*i);
        plot(t, v_RTN_1(:, i) - v_RTN_2(:, i), 'r', 'LineWidth', 1.5);
        ylabel(['v_{RTN} - ' titles{i} ' [km/s]']);
        if i == 1
            title('Error in RTN Velocity Components');
        end
        grid on;
    end
    
    xlabel('Time [s]');
    sgtitle(title_string);

    figure('Units', 'normalized', 'OuterPosition', [0 0 1 1]);  % Fullscreen
    
    titles = {'Radial (R)', 'Transverse (T)', 'Normal (N)'};
    
    for i = 1:3
        % Position subplot (left column)
        subplot(3, 2, 2*i - 1);
        plot(t, r_RTN_1(:, i), 'b', 'LineWidth', 1.5); hold on;
        plot(t, r_RTN_2(:, i), 'r--', 'LineWidth', 1.5);
        ylabel(['r_{RTN} - ' titles{i} ' [km]']);
        if i == 1
            title('RTN Position Components');
            legend(label1, label2);
        end
        grid on;
    
        % Velocity subplot (right column)
        subplot(3, 2, 2*i);
        plot(t, v_RTN_1(:, i), 'b', 'LineWidth', 1.5); hold on;
        plot(t, v_RTN_2(:, i), 'r--', 'LineWidth', 1.5);
        ylabel(['v_{RTN} - ' titles{i} ' [km/s]']);
        if i == 1
            title('RTN Velocity Components');
            legend(label1, label2);
        end
        grid on;
    end
    
    xlabel('Time [s]');
    sgtitle(title_string);
