function [new_D] = sequential_rejection_DAG(G, nodes_info)
% This is the main rejection algorithm to determine the rejected and
% accepted symmetry hypothesis

% input
% G   -- the graph
% nodes_info -- various nodes info (alpha, p values, etc.)


num_of_nodes = size(G.Nodes,1);
Alpha = 0.05;
leaves = zeros(num_of_nodes,1);
outDeg = outdegree(G);
leavesIdx = outDeg == 0;
leaves(leavesIdx) = 1;

% initializing
D = zeros(size(G.Nodes,1), 1);
counter = 1;
max_count = nodes_info(1).order; % chosen arbitrary for now

while counter<max_count
    % step 1: update the critical values
    weights = weight_update_for_sequential_rejection(G, D);
    num_leaves_not_in_D = leaves'*(1-D);

    for j=1:num_of_nodes
        nodes_info(j).alpha = (weights(j)*Alpha)/num_leaves_not_in_D;
    end

    % step 2: update the rejection set
    new_D = D;
    for j=1:num_of_nodes
        if and(D(j)==0,nodes_info(j).p_value<=nodes_info(j).alpha)
            new_D(j)=1;
        end
    end
    % reject all ancestor
    for j=1:num_of_nodes
        if new_D(j)==1
            ancestors = get_all_ancestors(G, j);
            new_D(ancestors)=1;
        end
    end

    if new_D==D
        % Done
        break;
    end
    % prepare for the next iteration
    counter = counter+1;
    D = new_D; 
end

end