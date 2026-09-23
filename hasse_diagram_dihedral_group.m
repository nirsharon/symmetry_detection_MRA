function [G, nodes_info] = hasse_diagram_dihedral_group(d)
    % Builds the Hasse diagram for Dihedral MRA (Collapsed Regime: odd d)
    % Based on Section 5.3.1, Case 1 of the paper.
    
    if mod(d, 2) == 0
        error('This function currently assumes d is odd (collapsed regime).');
    end
    
    % Get sorted divisors descending (e.g., [81, 27, 9, 3, 1])
    % Note: Replace 'find_divisors' with your actual function name if different
    divs = sort(find_divisors(d), 'descend'); 
    
    nodes_info = struct('order', {}, 'type', {}, 'alpha', {}, 'p_value', {});
    
    % Maps to easily look up node indices when building edges
    D_idx = containers.Map('KeyType', 'double', 'ValueType', 'double');
    C_idx = containers.Map('KeyType', 'double', 'ValueType', 'double');
    
    curr_idx = 1;
    
    % 1. Create all D_m nodes
    for i = 1:length(divs)
        m = divs(i);
        nodes_info(curr_idx).order = m;
        nodes_info(curr_idx).type = 'D';
        nodes_info(curr_idx).alpha = 0.05;
        nodes_info(curr_idx).p_value = 0;
        D_idx(m) = curr_idx;
        curr_idx = curr_idx + 1;
    end
    
    % 2. Create C_m nodes (excluding C_d and C_1 as per the paper)
    for i = 1:length(divs)
        m = divs(i);
        if m == d || m == 1
            continue;
        end
        nodes_info(curr_idx).order = m;
        nodes_info(curr_idx).type = 'C';
        nodes_info(curr_idx).alpha = 0.05;
        nodes_info(curr_idx).p_value = 0;
        C_idx(m) = curr_idx;
        curr_idx = curr_idx + 1;
    end
    
    edges = [];
    
    % 3. Build the directed edges (Implications)
    for i = 1:length(divs)
        m1 = divs(i);
        
        % D_{m1} -> D_{m2} (if m1/m2 is prime)
        for j = i+1:length(divs)
            m2 = divs(j);
            if mod(m1, m2) == 0 && isprime(m1 / m2)
                edges = [edges; D_idx(m1), D_idx(m2)];
            end
        end
        
        % C_{m1} -> C_{m2} (if m1/m2 is prime)
        if isKey(C_idx, m1)
            for j = i+1:length(divs)
                m2 = divs(j);
                if isKey(C_idx, m2) && mod(m1, m2) == 0 && isprime(m1 / m2)
                    edges = [edges; C_idx(m1), C_idx(m2)];
                end
            end
        end
        
        % D_{m1} -> C_{m1} (Dihedral implies Cyclic)
        if isKey(C_idx, m1)
            edges = [edges; D_idx(m1), C_idx(m1)];
        end
    end
    
    % Create the directed graph
    G = digraph(edges(:,1), edges(:,2));
end