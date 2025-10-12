function[] = plot_the_hasse_diagram(G, nodes_info)
% plotting the graph and the associated values on the nodes

to_save = 0;

num_subgroups = size(G.Nodes,1);
d = nodes_info(1).order;
figure;
h = plot(G, 'Layout', 'layered', 'Direction', 'down', ...
    'NodeLabel', {}, 'NodeColor', 'c', 'MarkerSize', 30, 'LineWidth', 1.5);
highlight(h, 'Edges', 1:numedges(G), 'EdgeColor', 'k', 'LineWidth', 1.2);

% Determine node colors: red for one entries in D, cyan otherwise (zero)
%node_colors = repmat([1 0 0], num_subgroups, 1);  % default cyan (RGB [0 1 1])
%zero_indices = find(D == 0);
%node_colors(zero_indices, :) = repmat([0 1 1], length(zero_indices), 1); % red

h.NodeColor = node_colors;

% Print subgroup labels inside nodes
for i = 1:num_subgroups
    x = h.XData(i);
    y = h.YData(i);
    text(x, y, sprintf('C_{%d}', nodes_info(i).order), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 10, 'FontWeight', 'bold');
end

if to_save
    folder_name = 'reject_process';
    if isfolder(folder_name)
        cd(folder_name);
    else
        mkdir(folder_name)
        cd(folder_name)
    end

    name_it = ['d_', num2str(d),'_Iteration_num_', num2str(iter_num) ] ;
    saveas(gcf, name_it ,'fig');
    saveas(gcf,name_it,'jpg');
    print('-depsc2',name_it);
    print('-depsc2',name_it);
    %   save('data_for_run');
    cd '../'
end
