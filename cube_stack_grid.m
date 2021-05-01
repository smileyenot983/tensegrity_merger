function [points] = cube_stack_grid(N,n_stacks,shift)

    n = n_stacks * (N^3);
    points = zeros(3,n);
    
    points(:,1:n/n_stacks) = cube_grid(N);
    
    step_size = n/n_stacks;
    for i=2:n_stacks
        
        points(:,(i-1)*step_size+1:i*step_size) = cube_grid(N);

        points(3,(i-1)*step_size+1:i*step_size) = points(3,(i-1)*step_size+1:i*step_size) + points(3,(i-1)*step_size) + ...
            (points(3,(i-1)*step_size) - points(3,(i-1)*step_size-1));
        

    end
    
    switch shift
        case 'x'
            points(1,:) = points(1,:) + 1 * max(points(1,:));
            disp("Shift in x");
        case 'y'
            points(2,:) = points(2,:) + 1 * max(points(2,:));
            disp("Shift in y");
        case 'z'
            points(3,:) = points(3,:) + 1 * max(points(3,:));
            disp("Shift in z");   
        otherwise
            disp("No shifting");     
    end
    
    
    
    
    
    
    