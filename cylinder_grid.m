function [points] = cylinder_grid(radius,n_points,shift)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%     In case shift = some coordinate ('x','y','z') => shifts corresponding
%     coordinates in a given shift direction by diameter(2*radius) of cylinder 

%     radius fixed
    
    
    angle = linspace(0,6.28,n_points/4);
    height = linspace(0,4,4);
    
    [theta, r, h] = meshgrid(angle, radius, height);

    [x,y,z] = pol2cart(theta,r,h);
    
    
    
    
    points(1,:) = x(:);
    points(2,:) = y(:);
    points(3,:) = z(:);
    
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