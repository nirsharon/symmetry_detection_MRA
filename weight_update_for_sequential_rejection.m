function [weights] = weight_update_for_sequential_rejection(G, D)
% This is where we update the weights for the rejection algorithm
%
% input
% node_info -- the structure of all nodes info
% D -- the current rejection set

num_of_nodes = size(G.Nodes,1);
weights = zeros(num_of_nodes,1);
edges_arr = table2array(G.Edges);

leaves = zeros(num_of_nodes,1);
outDeg = outdegree(G);
leavesIdx = outDeg == 0;
leaves(leavesIdx) = 1;

% initializing weights
for j=1:num_of_nodes
    if and(D(j)==0,leaves(j)==1)
        weights(j) = 1;
    else
        weights(j) = 0;
    end
end


% back loop from leaves to root
for j=num_of_nodes:-1:1
    if D(j)==0
        parents_j = edges_arr(edges_arr(:,2)==j,1);
        parents_num = numel(parents_j);
        parents_num_not_in_D = 0;
        for l=1:parents_num
            if D(parents_j(l))==0
            weights(parents_j(l)) = weights(parents_j(l)) + weights(j)/parents_num;
            parents_num_not_in_D = parents_num_not_in_D + 1;
            end
        end
        if parents_num>0
            weights(j) = weights(j) - (parents_num_not_in_D)*weights(j)/parents_num;
        end
    end
end

end