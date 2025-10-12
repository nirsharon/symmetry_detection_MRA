function [G, nodes_info] = DAG_dihedral_even(n)
% The constructor of the DAG diagram. 
%

    % Validate input
    if ~isscalar(n) || n < 1 || floor(n) ~= n
        error('Input must be a positive integer');
    end
    if mod(n,2)
        error('Input must be an even integer');
    end

    % Find divisors of n (each corresponds to a subgroup)
    divisors = sort(divisors_of_n(n), 'descend');
    num_subgroups = length(divisors);

    % Create nodes for each subgroup
    nodes_info = struct('order', {}, 'generator', {}, 'alpha', {}, 'p_value', {});
    for i = 1:num_subgroups
        m = divisors(i);
        g = mod(n/m, n); % generator of subgroup
        nodes_info(i).order = m;
        nodes_info(i).generator = g;
        nodes_info(i).is_dihedral = 1; 
        % adding cyclic node
        if and(m~=n,and(m~=1,m~=n))
            nodes_info(i+num_subgroups-2).order = m;
            nodes_info(i+num_subgroups-2).generator = g;
            nodes_info(i+num_subgroups-2).is_dihedral = 0; 
        end
    end

    % Build edges based on inclusion (minimal edges only)
    edges = [];
    for i = 1:num_subgroups
        for j = i+1:(num_subgroups)
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
                    if (j + num_subgroups-2)<=size(nodes_info,2)
                        edges = [edges; j (j + num_subgroups-2)]; % from Dm to Cm
                    end
                    %inside Cm, m<d
                    if and(i>2,(j + num_subgroups-2)<=size(nodes_info,2))
                        edges = [edges; (i + num_subgroups-2) (j + num_subgroups-2)];
                    end
                end
            end
        end
    end
    % form the tree graph
    G = digraph(edges(:,1), edges(:,2));

    % Extract edges
    E = table2array(G.Edges(:,1));  % get EndNodes

    % Remove duplicates
    [~, ia] = unique(sort(E,2), 'rows');

    % Rebuild graph with unique edges
    G = digraph(E(ia,1), E(ia,2));

end

%h = plot(G, 'Layout', 'layered', 'Direction', 'down', ...
%    'NodeLabel', {}, 'Marker', 'none', 'LineWidth', 1.5);
%highlight(h, 'Edges', 1:numedges(G), 'EdgeColor', 'k', 'LineWidth', 1.2);


