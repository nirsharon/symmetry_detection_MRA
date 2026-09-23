function [G, nodes_info] = hasse_diagram_dihedral_group_OLD(n)
% Build Hasse diagram (subgroup lattice) for D_{2n}, keeping your cyclic style.
% Labels follow your convention: D_n has order 2n; sub-dihedrals are D_m.

    if ~isscalar(n) || n < 3 || floor(n) ~= n
        error('Input must be an integer n >= 3 (so |D_{2n}| = 2n >= 6).');
    end

    divs = sort(divisors_of_n(n), 'descend');  % divisors of n

    nodes_info = struct('order', {}, 'generator', {}, 'alpha', {}, 'p_value', {});
    labels   = {};       % for NodeLabel
    kinds    = {};       % 'C','R*','D','G'
    params_m = [];       % m for C_m / D_{2m}

    % 1) Cyclic rotation subgroups C_m (m|n), descending m, SKIP m=1
    for m = divs
        if m == 1, continue; end               % <-- remove C_1
        nodes_info(end+1).order   = m;         % |C_m|=m
        nodes_info(end).generator = mod(n/m, n);
        nodes_info(end).alpha     = 0.05;
        nodes_info(end).p_value   = 0;
        labels{end+1}             = sprintf('C_{%d}', m);
        kinds{end+1}              = 'C';
        params_m(end+1)           = m;
    end

    % 2) Collapsed reflections node (all <r^k s>)
    nodes_info(end+1).order   = 2;
    nodes_info(end).generator = -1;            % placeholder
    nodes_info(end).alpha     = 0.05;
    nodes_info(end).p_value   = 0;
    labels{end+1}             = 'C_{2}^{(refl)}';  % <-- disambiguate
    kinds{end+1}              = 'R*';
    params_m(end+1)           = 0;

    % 3) Dihedral subgroups D_{2m} for m|n, m<n (descending m)
    for m = divs
        if m < n
            nodes_info(end+1).order   = 2*m;   % |D_{2m}|=2m
            nodes_info(end).generator = -m;    % placeholder
            nodes_info(end).alpha     = 0.05;
            nodes_info(end).p_value   = 0;
            labels{end+1}             = sprintf('D_{%d}', m);  % D_m (your style)
            kinds{end+1}              = 'D';
            params_m(end+1)           = m;
        end
    end

    % 4) Top group D_{2n}
    nodes_info(end+1).order   = 2*n;
    nodes_info(end).generator = 0;
    nodes_info(end).alpha     = 0.05;
    nodes_info(end).p_value   = 0;
    labels{end+1}             = sprintf('D_{%d}', n);          % D_n (your style)
    kinds{end+1}              = 'G';
    params_m(end+1)           = n;

    N = numel(nodes_info);

    % ------- edges (containment) then keep only covers ----------
    edges_all = [];

    % Cyclic chain: C_m -> C_d when d | m
    C_idx = find(strcmp(kinds,'C'));
    for ii = 1:numel(C_idx)
        i = C_idx(ii);  Mi = params_m(i);
        for jj = 1:numel(C_idx)
            j = C_idx(jj);  Dj = params_m(j);
            if Mi ~= Dj && mod(Mi, Dj) == 0
                edges_all(end+1,:) = [i, j]; %#ok<AGROW>
            end
        end
    end

    % D_{2m} includes C_d for d | m  → edge D_{2m} -> C_d
    D_idx = find(strcmp(kinds,'D'));
    for ii = 1:numel(D_idx)
        i = D_idx(ii);  Mi = params_m(i);
        for jj = 1:numel(C_idx)
            j = C_idx(jj);  Dj = params_m(j);
            if mod(Mi, Dj) == 0
                edges_all(end+1,:) = [i, j]; %#ok<AGROW>
            end
        end
    end

    % Dihedral chain: D_{2m} -> D_{2d} when d | m, d < m
    for ii = 1:numel(D_idx)
        i = D_idx(ii);  Mi = params_m(i);
        for jj = 1:numel(D_idx)
            j = D_idx(jj);  Dj = params_m(j);
            if Mi ~= Dj && mod(Mi, Dj) == 0
                edges_all(end+1,:) = [i, j]; %#ok<AGROW>
            end
        end
    end

    % Reflections contained in every D_{2m} and in the top D_{2n}
    R_idx = find(strcmp(kinds,'R*'));
    for ii = 1:numel(D_idx)
        i = D_idx(ii);
        edges_all(end+1,:) = [i, R_idx]; %#ok<AGROW>
    end
    G_idx = find(strcmp(kinds,'G'));
    edges_all(end+1,:) = [G_idx, R_idx];

    % Top includes all D_{2m} and all C_m
    for ii = 1:numel(D_idx), edges_all(end+1,:) = [G_idx, D_idx(ii)]; end
    for ii = 1:numel(C_idx), edges_all(end+1,:) = [G_idx, C_idx(ii)]; end

    % ---- keep only cover edges (minimal) ----
    edges_all = unique(edges_all, 'rows');
    A = false(N,N);
    for e = 1:size(edges_all,1), A(edges_all(e,1), edges_all(e,2)) = true; end
    keep = true(size(edges_all,1),1);
    for e = 1:size(edges_all,1)
        i = edges_all(e,1); j = edges_all(e,2);
        for k = 1:N
            if k~=i && k~=j && A(i,k) && A(k,j), keep(e)=false; break; end
        end
    end
    edges_cover = edges_all(keep,:);

    G = digraph(edges_cover(:,1), edges_cover(:,2), [], N);

    % Optional labels for plotting
    try, G.Nodes.Label = labels(:); end
end

function divs = divisors_of_n(n)
    divs = [];
    for k = 1:n
        if mod(n,k) == 0, divs = [divs k]; end %#ok<AGROW>
    end
end
