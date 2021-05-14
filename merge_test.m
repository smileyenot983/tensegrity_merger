clc;
close all;
clear all;
cvx_clear;
cvx_solver Gurobi_2;


%  in case you don't have any pregenerated structures
sol = run_experiment(2,10,true,"cylinder",0,0,false,false,...
                        false,0);

save("random_structure1","sol");
                    
sol = run_experiment(2,10,true,"cylinder",0,0,false,false,...
                        false,0);
                    
save("random_structure2","sol");


%  load structures
path1 = "random_structure1.mat";
path2 = "random_structure2.mat";

sol1 = load("random_structure1","sol");
sol2 = load("random_structure2","sol");


% merge
shift_direction = "x";
shift_distance = 2.5;

seed = 1;
merge_result = merge(seed,path1,path2,shift_direction,shift_distance);










