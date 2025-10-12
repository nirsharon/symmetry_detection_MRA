function [G, nodes_info] = hasse_diagram_cyclic_group(n)
% The constructor of the hasse diagram. 
% initial alpha value is hardcoded to be 0.05
% initial p value is set to 0

    % Validate input
    if ~isscalar(n) || n < 1 || floor(n) ~= n
        error('Input must be a positive integer');
    end

    % Find divisors of n (each corresponds to a subgroup)
    divisors = sort(divisors_of_n(n), 'descend');
    num_subgroups = length(divisors)-1;

    % Create nodes for each subgroup
    nodes_info = struct('order', {}, 'generator', {}, 'alpha', {}, 'p_value', {});
    for i = 1:num_subgroups
        d = divisors(i);
        g = mod(n/d, n); % generator of subgroup
        nodes_info(i).order = d;
        nodes_info(i).generator = g;
        nodes_info(i).alpha = 0.05; %rand();    % Replace ?
        nodes_info(i).p_value = 0; %rand();     % Replace ?
    end

    % Build edges based on inclusion (minimal edges only)
    edges = [];
    for i = 1:num_subgroups
        for j = i+1:num_subgroups
            if mod(divisors(i), divisors(j)) == 0
                is_minimal = true;
                for k = i+1:j-1
                    if mod(divisors(i), divisors(k)) == 0 && mod(divisors(k), divisors(j)) == 0
                        is_minimal = false;
                        break;
                    end
                end
                if is_minimal
                    edges = [edges; i j]; % Directed edge from larger to smaller subgroup
                end
            end
        end
    end
    % form the tree graph
    G = digraph(edges(:,1), edges(:,2));
end
