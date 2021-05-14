function sol_merg = merge(seed,path1,path2,shift,shift_dist)
%MERGE Summary of this function goes here
%   path should contain structure with following parameters:
% R,C,points 

foldername = strcat("Merge results. Seed:",string(seed));
mkdir(foldername);

sol1 = load(path1);
sol2 = load(path2);

sol1 = sol1.sol;
sol2 = sol2.sol;

switch shift
    case 'x'
        sol2.points(1,:) = sol2.points(1,:) + shift_dist;
        disp("Shift in x");
    case 'y'
        sol2.points(2,:) = sol2.points(2,:) + + shift_dist;
        disp("Shift in y");
    case 'z'
        sol2.points(3,:) = sol2.points(3,:) + + shift_dist;
        disp("Shift in z");        
    otherwise
        disp("No shifting");
end


%% Combining together 
C_1 = sol1.C;
C_2 = sol2.C;

C_bar = zeros(size(C_1,1)+size(C_2,1),size(C_1,2)+size(C_2,2));

C_bar(1:size(C_1,1),1:size(C_1,2)) = C_1;
C_bar(size(C_1,1)+1:end,size(C_1,2)+1:end) = C_2;


% concatenating strut connectivity into diagonal matrix
R_1 = sol1.R;
R_2 = sol2.R;

R_bar = zeros(size(R_1,1)+size(R_2,1),size(R_1,2)+size(R_2,2));
                
R_bar(1:size(R_1,1),1:size(R_1,2)) = R_1;
R_bar(size(R_1,1)+1:end,size(R_1,2)+1:end) = R_2;


% concatenating node coordinates horizontally

p_1 = sol1.points;
p_2 = sol2.points;

p_bar = horzcat(p_1,p_2);



%% Structure before merging

% big M relaxation variable:
big_M = 20;

% number of fixed nodes:
m=3;

% total number of nodes:
new_n = size(R_bar,1);

% gravity
external_force = repmat([0;0;-9.8], [1, new_n]);    


Dir = zeros(new_n,new_n,3);
for i =1:new_n
    for j =1:new_n
        
        Dir(i,j,:) = reshape(  (p_bar(:,i) - p_bar(:,j)), [1, 1, 3]  );
        
    end
end

Dir_t = permute(Dir, [3, 2, 1]);

separate_sol.points = p_bar;
separate_sol.R = R_bar;
separate_sol.C = C_bar;

dataname = strcat("unmerged_",".mat");


filename = strcat(foldername,"/",dataname);

visualize_2(separate_sol,3,filename,"separate")


%% OPtimization

strut_constraint = false;
cable_constraint = false;
projection_constraint = false;

cvx_begin

variable delta_R(new_n, new_n) binary
variable delta_C(new_n, new_n) binary

variable f(new_n, new_n)
variable g(3, m)

minimize(sum(delta_C(:)) + sum(delta_R(:)))

C = C_bar + delta_C;
R = R_bar + delta_R;
points = p_bar;

subject to


    C == C';
    R == R';

    sum(C,1) >= 3*ones(1,new_n);
    
    diag(C) == zeros(new_n,1);
    diag(R) == zeros(new_n,1);
    
    P = C + R;
    P(:) <= ones(new_n*new_n,1);
    
    
    for i=1:new_n
        for j=1:new_n
            f(i,j) <= big_M * C(i,j); 
            -f(i,j) <= big_M * R(i,j);
        end 
    end
    
    g_ext = [g,zeros(3,new_n-m)];
    
    
    for i=1:new_n
        Dir_t(:,:,i) * f(i,:)' == external_force(:,i) + g_ext(:,i);
       
    end
    
    for i=1:new_n
        Dir_t(:,:,i) * f(i,:)' == external_force(:,i);
       
    end


%     slice part, which is responsible for connection between 2 separate
%     structures
%     C_slice = delta_C(1:n,n+1:new_n);
%     sum(C_slice(:)) >= min_C; 
%     
% 
%     R_slice = delta_R(1:n,n+1:new_n);
%     sum(R_slice(:)) >= min_R; 
    
    
    if strut_constraint==true
        for i=1:new_n        
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

cvx_end


sol_merg.points = p_bar;
sol_merg.R = R_bar+delta_R;
sol_merg.C = C_bar+delta_C;
sol_merg.delta_R = delta_R;
sol_merg.delta_C = delta_C;

dataname = strcat("merged_",".mat");

filename = strcat(foldername,"/",dataname);

% save(filename,"sol1","sol2","sol_merg");

visualize_2(sol_merg,3,filename,"merged");





end

