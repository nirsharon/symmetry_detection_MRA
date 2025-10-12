% script name: "plot_hypothesis_DAG"

close all;
clear;

% basic parameters
d = 72;

to_save = 0;

[G, nodes_info] = hasse_diagram_cyclic_group(d);

num_subgroups = size(G.Nodes,1);
    
figure;
    
% Plot the graph without nodes
h = plot(G, 'Layout', 'layered', 'Direction', 'down', ...
    'NodeLabel', {}, 'Marker', 'none', 'LineWidth', 1.5);
highlight(h, 'Edges', 1:numedges(G), 'EdgeColor', 'k', 'LineWidth', 1.2);

% Add transparent circles for nodes
hold on;
scatter(h.XData, h.YData, 600, 'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', 'w', 'LineWidth', 1.5);
hold off;

% Determine node colors: red for one entries in D, cyan otherwise (zero)
%node_colors = repmat([1 0 0], num_subgroups, 1);  % default cyan (RGB [0 1 1])
%zero_indices = find(D == 0);
%node_colors(zero_indices, :) = repmat([0 1 1], length(zero_indices), 1); % red
% h.NodeColor = node_colors;

% Print subgroup labels inside nodes
for i = 1:num_subgroups
    x = h.XData(i);
    y = h.YData(i);
    text(x, y, sprintf('C_{%d}', nodes_info(i).order), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 10, 'FontWeight', 'bold');
end


if to_save
    folder_name = 'DAG';
    if isfolder(folder_name)
        cd(folder_name);
    else
        mkdir(folder_name)
        cd(folder_name)
    end

    name_it = ['hasse_d_', num2str(d)] ;
    saveas(gcf, name_it ,'fig');
    saveas(gcf,name_it,'jpg');
    exportgraphics(gcf, [name_it,'.pdf'], 'ContentType', 'vector');
%   save('data_for_run');
    cd '../'
end
