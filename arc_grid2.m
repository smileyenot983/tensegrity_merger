function [points] = arc_grid2(n_rad,n_wid,n_ang,shift);
% 

% n_rad=2;
% n_wid=2;
% n_ang = 10;

    points = zeros(3,n_rad*n_wid*n_ang);
    angle = linspace(0,pi,n_ang);

    min_rad=1;
    last_index=1;
    for i=1:n_rad
        x = (i-1+min_rad)*cos(angle); 
        y = (i-1+min_rad)*sin(angle);

        for j=1:n_ang

            for k = 1:n_wid
                points(:,last_index) = [x(j);y(j);k];
                last_index =last_index+1;
            end

        end

    end


    switch shift
        case 'x'
            points(1,:) = points(1,:) + 2 * max(points(1,:));
            disp("Shift in x");
        case 'y'
            points(2,:) = points(2,:) + 2 * max(points(2,:));
            disp("Shift in y");
        case 'z'
            points(3,:) = points(3,:) + 2 * max(points(3,:));
            disp("Shift in z");        
        otherwise
            disp("No shifting");
    end






end

