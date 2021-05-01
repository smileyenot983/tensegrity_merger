function [points] = sphere_grid(radius,n_points,shift)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    % Define the radius of the shell of points, the surface of the sphere.
%     radius = 10; 
    % Define the number of points to place on the surface of the sphere.
    numPoints = n_points;	% Use a large value.
    % Get a 3-by-numPoints list of (x, y, z) coordinates.
    r = randn(3, numPoints);
    % At this point the points can be anywhere in 3D space,
    % not just on the surface of a sphere, a shell of constant radius.
    % Divide by the distance of each point from the origin.
    % This will place each point out at a definite radius of 1, not scattered throughout the volume.
    r = bsxfun(@rdivide, r, sqrt(sum(r.^2,1)));
    % Now multiply by radius to make the shell out 
    % at a specified radius instead of a radius of exactly 1
    points = radius * r;
    % Extract the x, y, and z coordinates from the array.
%     x = r(1,:); % Extract x from row #1.
%     y = r(2,:); % Extract y from row #2.
%     z = r(3,:); % Extract z from row #3.
%     % Display the shell of points
%     scatter3(x, y, z);
%     axis square; % Make sure the aspect ratio is maintained as it's displayed and rotated.
%     xlabel('X', 'FontSize', 20);
%     ylabel('Y', 'FontSize', 20);
%     zlabel('Z', 'FontSize', 20);
    % Enlarge figure to full screen.
    % set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    % msgbox('Now use the circular arrow icon on the toolbar to rotate the sphere.');

    
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

