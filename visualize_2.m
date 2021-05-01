function [] = visualize_solution(solution,n_fixed,filename,plot_title)

    FontName = 'Tibetan Machine Uni';

%VISUALIZE_GENERATION Summary of this function goes here
%   Detailed explanation goes here
    figure1 = figure;
    points = solution.points;
    C = solution.C;
    R = solution.R;
    n = size(points,2);
    
    
    for i=1:n
        P = points(:,i);
        
        text(P(1), P(2), P(3), num2str(i) ,'HorizontalAlignment','left','FontSize',16,'FontName',FontName);
        for j=1:n
            if C(i,j)==1
                P1 = points(:,i);
                P2 = points(:,j);
                plot3([P1(1), P2(1)], [P1(2), P2(2)], [P1(3), P2(3)], ...
                    'LineWidth', 1, 'Color', [0,0,0]); hold on;
            
            end
            
            if R(i,j)==1
                P1 = points(:,i);
                P2 = points(:,j);
%                 plot3([P1(1), P2(1)], [P1(2), P2(2)], [P1(3), P2(3)], ...
%                     'LineWidth', 2, 'Color', [100,0,0]/255); hold on;
                vis_Cylinder(P1, P2, 0.05, ...
                    'EdgeAlpha', 0, 'FaceAlpha', 1, 'FaceColor', [0.9 0.2 0.3], 'SpecularStrength', 0.2); hold on;
                
            end

        end
        
    end
    
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    
    for i = 1:size(points, 2)
        if i > n_fixed
            PointColor = [0.2 1 0.5];
        else
            PointColor = [0 0.1 1];
        end
    
    vis_Sphere(points(:, i)', 0.15, ...
        'EdgeAlpha', 0, 'FaceAlpha', 1, 'FaceColor', PointColor, 'SpecularStrength', 0.2)
    end
    
    
%     plot3(points(1, :)', points(2, :)', points(3, :)', 'o', 'MarkerFaceColor', 'r')
%     plot3(points(1, 1:n_fixed)', points(2, 1:n_fixed)', points(3, 1:n_fixed)', 'o', 'MarkerFaceColor', 'g')

%     title(plot_title);
%     axis equal
%     
%     saveas(figure1,strcat(filename,'.png'));

    grid on; grid minor;
ax = gca;
ax.GridAlpha = 0.6;
ax.LineWidth = 0.5;
ax.MinorGridLineStyle = '-';
ax.MinorGridAlpha = 0.2;
ax.FontName = FontName;
% ax.FontName = 'Tibetan Machine Uni';
ax.FontSize = 14;
xlabel_handle = xlabel('$$x$$, m', 'Interpreter', 'latex');
ylabel_handle = ylabel('$$y$$, m', 'Interpreter', 'latex');
zlabel_handle = zlabel('$$z$$, m', 'Interpreter', 'latex');

axis equal
camlight('right')

    drawnow;
%     mkdir(['./', convertStringsToChars(path), '_paper']);
%     filename = ['./', convertStringsToChars(path), '_paper', '/',solution_num, '_', convertStringsToChars(plot_title)];
%     saveas(gcf, [filename, '.png'], 'png');
    
    saveas(gcf,strcat(filename,'.png'));
    
%     savefig([filename, '.fig']);
%     savefig(filename)
    
end

