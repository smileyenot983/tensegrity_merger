function [] = visualize_solution(solution,n_fixed,filename,plot_title)
%VISUALIZE_GENERATION Summary of this function goes here
%   Detailed explanation goes here
    figure1 = figure;
    points = solution.points;
    C = solution.C;
    R = solution.R;
    n = size(points,2);
    for i=1:n
        P = points(:,i);
        
        text(P(1), P(2), P(3), num2str(i) ,'HorizontalAlignment','left','FontSize',16);
        for j=1:n
            if C(i,j)==1
                P1 = points(:,i);
                P2 = points(:,j);
                plot3([P1(1), P2(1)], [P1(2), P2(2)], [P1(3), P2(3)], ...
                    'LineWidth', 1, 'Color', [189,183,1]/255); hold on;
            
            end
            
            if R(i,j)==1
                P1 = points(:,i);
                P2 = points(:,j);
                plot3([P1(1), P2(1)], [P1(2), P2(2)], [P1(3), P2(3)], ...
                    'LineWidth', 2, 'Color', [100,0,0]/255); hold on;
                
            end

        end
        
    end
    
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    
    plot3(points(1, :)', points(2, :)', points(3, :)', 'o', 'MarkerFaceColor', 'r')
    plot3(points(1, 1:n_fixed)', points(2, 1:n_fixed)', points(3, 1:n_fixed)', 'o', 'MarkerFaceColor', 'g')

    title(plot_title);
    axis equal
    
    saveas(figure1,strcat(filename,'.png'));
    
end

