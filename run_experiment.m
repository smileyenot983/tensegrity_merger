function solution = run_experiment(exp_num,n_nodes,random_nodes,grid_structure,shift,shift_dist,strut_constraint,cable_constraint,...
    projection_constraint,projection_axis)
%RUN_EXPERIMENT Summary of this function goes here
%   Detailed explanation goes here:

%   exp_num  -  seed [int]
%   n_nodes  -  number of nodes for structure [int]
%   random_nodes  -  whether to choose n_nodes randomly or not [bool]
%   grid_structure  -  shape of points from which we sample [str]
%   shift  -  shift in given direction(axis) [char]

%   strut_constraint  -  whether to use constraint on strut max length [bool]
%   strut_max_length  -  maximum length of strut [float]

%   cable_constraint  -  whether to use constraint on cable max length [bool]
%   cable_max_length  -  maximum length of cable [float]

%   projection_constraint  -  whether to use constraint on maximum projection [bool]
%   projection_axis  -  axis of max projection constraint [str]
%   projection_max_length  -  maximum length of a projection

%% GRID

rng(exp_num);
switch grid_structure
    case "cube stack"
        N=3;
        stacks = 2;
        points = cube_stack_grid(N,stacks,shift,shift_dist);
        
        strut_max_length = 4; 
        cable_max_length = 4; 
        projection_max_length = 4;
    case "sphere"
        n = 100;
        radius = 3;
        points = sphere_grid(radius,n,shift);
        
        strut_max_length = 4; 
        cable_max_length = 4; 
        projection_max_length = 4;
    case "cylinder"
        n = 100;
        radius = 2;
        points = cylinder_grid(radius,n,shift);
        
        strut_max_length = 4;
        cable_max_length = 4;
        projection_max_length = 4;
    case "arc"        
        n_rad=2;
        n_wid=3;
        n_ang = 10;
        points = arc_grid2(n_rad,n_wid,n_ang);
        
        strut_max_length = 6; 
        cable_max_length = 6; 
        projection_max_length = 6;
    otherwise
        points = randn(3,n);
        
end

% plotting grid structure
% scatter3(points(1,:),points(2,:),points(3,:));


%%

% after building grid of some shape we choose n points in this grid, to
% generate tensegrity
n = size(points,2);
m = 3;
% n_nodes= 20;
available_slots = 1:n;



if random_nodes == true
    random_permutation = randperm(numel(available_slots));
else
    load("cylinder_permutation");
    
end


chosen_slots = random_permutation(1:n_nodes);


chosen_points = points(:,chosen_slots);
n = n_nodes;
points = chosen_points;

Dir = zeros(n, n, 3);
for i = 1:n
    for j = 1:n
        
        Dir(i, j, :) = reshape( (points(:, i) - points(:, j)), [1, 1, 3] );
        
    end
end
Dir_t = permute(Dir, [3, 2, 1]);



% figure('Color', 'w', 'Name', 'arrow array')
% for i = 1:n
%     for j = 1:n
%         
%         arrow3d(...
%             'start', [i; j; 0], ...
%             'stop',  [i; j; 0] + reshape(Dir(i, j, :), [3, 1, 1]), ...
%             'ang', 15); 
%         hold on;
%         
%     end
% end
% 
% axis equal;
   
%% Checking distances between points

% norms = zeros(n);
% for i=1:n
%     
%     for j=1:n
%         norms(i,j) = norm(points(:,i)-points(:,j));
%         
%     end
%     
% end
%%

% constraint on max length of strut

external_force = repmat([0;0;-9.8], [1, n]);    

big_M = 10;

cvx_begin

variable R(n, n) binary
variable C(n, n) binary

variable f(n, n)
variable g(3, m)

variable c_1(n) binary
variable c_2(n) binary

variable c_3 binary
variable c_4 binary

minimize( sum(C(:)) )
subject to

    C == C';
    R == R';
    
    diag(R) == zeros(n, 1);
    diag(C) == zeros(n, 1);
    
%     to make at least 3 cable connections from each node
    sum(C,1) == 3*ones(1,size(C,1))
    
    
    sum(R, 1) == ones(1, n);
    sum(R, 2) == ones(n, 1);
    
    P = C+R;
    P(:) <= ones(n*n, 1);
    
%     !!!increasing number of cables!!!
%     sum(P(:)) >= 45;
    
   
    for i = 1:n
        for j = 1:n
           f(i, j) <= big_M * C(i, j);
          -f(i, j) <= big_M * R(i, j);
        end
    end
    
    g_ext = [g, zeros(3, n-m)];
    
    for i = 1:n 
        Dir_t(:, :, i) * f(i, :)' == external_force(:, i) + g_ext(:, i);
    end
    
    
    if strut_constraint==true
        for i=1:n        
            norm(sum(  repmat(R(i,:),3,1).*points,2  )-points(:,i)) <=strut_max_length;
        end
    end
    
    
    if cable_constraint==true
        for i=1:n
            for j = 1:n
                norm((points(:,i)-points(:,j)).*repmat(C(i,j),3,1)) <= cable_max_length;

            end

        end
    end

    x = [1;0;0];
    y = [0;1;0];
    z = [0;0;1];
    % constraint for restricting in one direction
    
    if projection_constraint == true
    
        if projection_axis=='x'
            dir=x;
        elseif projection_axis =='y'
            dir=y;
        elseif projection_axis =='z'
            dir=z;
        elseif projection_axis =='xy'
            dir=[x,y];
        else
            dir=[0;0;0];
        end


    %     l_0 = 2.5;
        for j = 1:size(dir,2)
            for i=1:n
                norm(  (sum(repmat(R(i,:),3,1).*points,2) -points(:,i))'*dir(:,j) ) <= projection_max_length;
            end
        end
    
    end
       
    
%     for i=1:n
%         norm(  points(:,i) - sum(repmat(R(i,:),3,1).*points,2)  ) <= l_0 + y'*(points(:,i) - sum(repmat(R(i,:),3,1).*points,2)) + big_M*(1-c_1(i));
%         norm(  points(:,i) - sum(repmat(R(i,:),3,1).*points,2)  ) <= l_0 - y'*(points(:,i) - sum(repmat(R(i,:),3,1).*points,2)) + big_M*(1-c_2(i));
%         c_1(i)+c_2(i) == ones(n,1);
%     end


cvx_end

C = full(C);
R = full(R);
f = full(f);

disp("sum: ");
disp(sum(P(:)));

%% Checkin directional norms
% x = [1;0;0];
% y = [0;1;0];
% z = [0;0;1];
% 
% z_norms = zeros(i); 
% for i=1:n
% %     z_norms(i) = norm(sum(  repmat(R(i,:),3,1).*points,2)'*z );
%     z_norms(i) = norm(  (sum( repmat(R(i,:),3,1).*points,2) - points(:,i))'*z     );
% 
% end

%% Plotting lengths of all connections:
plot_results=false;
if plot_results==true
    fprintf("_____________________________________________\n")
    disp("Strut lengths: ")
    for i=1:n
        for j=1:n
            if R(i,j)==1
                fprintf("Connection %d : %d  length:%f \n",i,j,norm(points(:,i)-points(:,j)));
    %             norm(points(:,i)-points(:,j))

            end

        end
        norm(sum(repmat(R(i,:),3,1).*points,2)-points(:,i));

    end
    fprintf("\n");

    %% Plotting lengths of all connections:

    fprintf("_____________________________________________\n")
    disp("Cable lengths: ")
    for i=1:n
        for j=1:n
            if C(i,j)==1
                fprintf("Connection %d : %d  length:%f \n",i,j,norm(points(:,i)-points(:,j)));
    %             norm(points(:,i)-points(:,j))

            end

        end
        norm(sum(repmat(R(i,:),3,1).*points,2)-points(:,i));

    end
    fprintf("\n");

    %% Plotting norms of all projections
    disp("Projection norms: ")
    for i=1:n
        for j=1:n
            if R(i,j)==1
                proj_x = norm(x'*(points(:,i)-points(:,j)));
                proj_y = norm(y'*(points(:,i)-points(:,j)));
                proj_z = norm(z'*(points(:,i)-points(:,j)));
                fprintf("Connection %d : %d, proj on x:%f, proj on y:%f, proj on z:%f \n",i,j,proj_x,proj_y,proj_z);
    %             norm(points(:,i)-points(:,j))

            end

        end
    %     norm(sum(repmat(R(i,:),3,1).*points,2)-points(:,i));

    end

    fprintf("_____________________________________________\n")
end
    %%
% Checking if all nodes are connected
% Using depth-first search(DFS)

check_connectivity=false;

if check_connectivity
    current_graph = graph(C + R);
    connection_check = dfsearch(current_graph,1);
    if size(connection_check,1) ~= n
        warning("Not all nodes are connected");
    end
end

%% Saving into a file

solution.points = points;
solution.C = C; 
solution.R = R;

% foldername = "last_tests";
% dataname = "nice_shape_1.mat";
% const = "unconstrained";
% if length_constraint==true
%     const = "constrained";
% end

% dataname = strcat(grid_structure,"_",exp_num,"_","seed_",string(seed),"_",const,".mat");
% dataname = strcat(grid_structure,"_",exp_num,"_",const,".mat");
% dataname = strcat(grid_structure,"_",string(exp_num),"_","c_constr:",string(cable_constraint),"_","s_constr:",string(strut_constraint),...
%     "_","proj_constr:",string(projection_constraint),".mat");
% mkdir(foldername);

% filename = "3_cables/nice_shape_2.mat";
% filename = strcat(foldername,"/",dataname);
% 
% save(filename,"solution");

% load(filename);
% points = solution.points;
% C = solution.C;
% R = solution.R;

%% 

% visualize = false;
% if visualize
%     visualize_solution(solution,m,filename);
% end
end

