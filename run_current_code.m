% script name: "run_current_code"

close all
clear

% main parameters
d = 20;
m = 5;




% construct the graph
[G, nodes_info] = hasse_diagram_cyclic_group(d);

% initial p_values
new_nodes_info = calc_raw_p_vals(nodes_info, synthetic_power_spec, 0.05, 0.5, 10000);

% plot 
plot_the_hasse_diagram(G, new_nodes_info)


[new_D] = sequential_rejection_DAG(G, nodes_info)


%===== old testing of weight_update_for_sequential_rejection
% D = randi([0, 1], size(G.Edges,1), 1);
% D = zeros(size(G.Edges,1), 1);
%[weights] = weight_update_for_sequential_rejection(G, D);
%=====
