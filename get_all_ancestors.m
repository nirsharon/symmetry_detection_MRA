function ancestor_list = get_all_ancestors(G, node_idx)
    % Initialize ancestor list
    ancestor_list = [];

    % Direct parents of the current node
    direct_parents = predecessors(G, node_idx);
    
    % For each direct parent, recursively find their ancestors
    for parent = direct_parents'
        ancestor_list(end+1) = parent; %#ok<AGROW>
        ancestor_list = [ancestor_list, get_all_ancestors(G, parent)]; %#ok<AGROW>
    end
    
    % Remove duplicates (if any)
    ancestor_list = unique(ancestor_list);
end