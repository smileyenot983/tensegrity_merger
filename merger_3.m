clc;
close all;
clear all;
cvx_clear;
cvx_solver Gurobi_2;

% THIS CODE IS WORKING
% 1.you can change hyperparameters of separate structure generation on line 17
% 2.you can change hyperparameters of merging on line 113


% in case separate structures are being generated, in order to avoid generating 2 
% identical structures for second structure seed = seed+1 
seed = 7;
foldername = strcat("Merge results. Seed:",string(seed));
mkdir(foldername);



%% PART 1: Obtaining 2 separate structures
% loading existing structures
load_structures = false;
% if existing structures are in separate files
load_separate = false;

path1 = "Merge results. Seed:12/cylinder_merged_c_constr:false_s_constr:true_proj_constr:false_ax:0";
path2 = "Merge results. Seed:12/cylinder_merged_c_constr:false_s_constr:true_proj_constr:false_ax:0";

% in this case there should be paths for both structures
if load_structures == true && load_separate==true
    sol1 = load(path1);
    sol2 = load(path2);

% in this case there should be a struct with fields: sol1, sol2
elseif load_structures == true && load_separate==false
    sol = load(path1);
    sol1 = sol.sol1;
    sol2 = sol.sol2;
else
    
    %______________PARAMETERS FOR GENERATING 2 SEPARATE STRUCTURES_________________
    grid_structure1 = "sphere";
    grid_structure2 = "cylinder";
    % number of nodes:
    n = 10;


    % shift direction(shift distance for cylinder = 2*radius)
    shift1 = 0;
    shift2 = 'x';


    % constraints [strut max length, cable max length, projection max length]
    constraints_set1 = [false,false,false];
    projection_axis1 = 0;

    constraints_set2 = [false,false,false];
    projection_axis2 = 0;
    
    %______________________________________________________________________________
    
    strut_constraint = constraints_set1(1);
    cable_constraint = constraints_set1(2);
    projection_constraint = constraints_set1(3);
    projection_axis = projection_axis1;
    
    sol1 = run_experiment(seed,n,true,grid_structure1,shift1,strut_constraint,cable_constraint,...
                    projection_constraint,projection_axis);


    strut_constraint = constraints_set2(1);
    cable_constraint = constraints_set2(2);
    projection_constraint = constraints_set2(3);
    projection_axis = projection_axis2;

    % second structure is shifted      
    sol2 = run_experiment(seed+1,n,true,grid_structure2,shift2,strut_constraint,cable_constraint,...
                        projection_constraint,projection_axis);

                
end


%% Preparations before merging             
% concatenating cable connectivity into diagonal matrix
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


%% PART 2: MERGING 2 SEPARATE STRUCTURES INTO SINGLE STRUCTURE


%_______________VARIABLES FOR OPTIMIZATION_________________

% big M relaxation variable:
big_M = 20;

% number of fixed nodes:
m=3;

% minimal number of cables between 2 separate structures
min_C = 1;
% minimal number of struts between 2 separate structures
min_R = 1;


% constraints for merging:
constraints_set_merg = [true,false,false];
projection_axis_merg = 0;

strut_max_length = 6;
cable_max_length = 6;
projection_max_length = 6;

%___________________________________________________________


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


dataname = strcat("unmerged_","c_constr:",string(cable_constraint),"_","s_constr:",string(strut_constraint),...
    "_","proj_constr:",string(projection_constraint),"_ax:",string(projection_axis),".mat");

filename = strcat(foldername,"/",dataname);

% visualize_solution(separate_sol,3,filename,"separate")
% visualize_solution(separate_sol,3,filename,"separate")

visualize_2(separate_sol,3,filename,"separate")
%%
% IT'S TIME TO OPTIMIZE THINGS

strut_constraint = constraints_set_merg(1);
cable_constraint = constraints_set_merg(2);
projection_constraint = constraints_set_merg(3);


% C = C_bar + delta_C
% R = R_bar + delta_R

% min || delta_C + delta_R ||

cvx_begin

variable delta_R1(new_n, new_n) binary
variable delta_C1(new_n, new_n) binary

variable f(new_n, new_n)
variable g(3, m)

minimize(sum(delta_C1(:)) + sum(delta_R1(:)))

C = C_bar + delta_C1 ;
R = R_bar + delta_R1 ;
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
sol_merg.R = R_bar + delta_R1 ;
sol_merg.C = C_bar + delta_C1 ;


%% SAVING RESULTS

dataname = strcat("merged_","c_constr:",string(cable_constraint),"_","s_constr:",string(strut_constraint),...
    "_","proj_constr:",string(projection_constraint),"_ax:",string(projection_axis),".mat");

filename = strcat(foldername,"/",dataname);

save(filename,"sol1","sol2","sol_merg");


% add image saving

% visualize_solution(sol_merg,3,filename,"merged");
% visualize_solution(sol_merg,3,filename,"merged");
visualize_2(sol_merg,3,filename,"merged");


%% Testing reducing

cvx_clear;

cvx_begin

variable delta_R2(new_n, new_n) binary
variable delta_C2(new_n, new_n) binary

variable f(new_n, new_n)
variable g(3, m)

maximize(sum(delta_C2(:)) + sum(delta_R2(:)))

C = C_bar - delta_C2;
R = R_bar - delta_R2;
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
    
    
    

cvx_end


sol_merg2.points = p_bar;

sol_merg2.R = sol_merg.R - delta_R2;
sol_merg2.C = sol_merg.C - delta_C2;




%% SAVING RESULTS

% dataname = strcat("merged_","c_constr:",string(cable_constraint),"_","s_constr:",string(strut_constraint),...
%     "_","proj_constr:",string(projection_constraint),"_ax:",string(projection_axis),".mat");
% 
% filename = strcat(foldername,"/",dataname);
% 
% save(filename,"sol1","sol2","sol_merg");


% add image saving

% visualize_solution(sol_merg,3,filename,"merged");
% visualize_solution(sol_merg,3,filename,"merged");
visualize_solution(sol_merg2,3,filename,"merged");
