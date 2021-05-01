function [points] = cube_grid(N)


%     N - number of nodes in each direction
    
    
    % CUBIC GRID
    % number of nodes in each direction
%     N = 3;
    % total number of nodes
    n = N^3;
    points = zeros(3,n);
    x_0 = 0;
    x_f = 2;
%     step_size = x_f/N;

    % number of fixed nodes
%     m = 3;

%     x=0;y=0;z=0;

    
    pos = 1;

    coords = linspace(x_0,x_f,N);
    for i=1:N
        for j=1:N
            for k=1:N

                x=coords(i);y=coords(j);z=coords(k);

                points(:,pos) = [x,y,z];
                pos=pos+1;
            end
        end
    end
end

