% script name: "plot_hypothesis_DAG_dihedral"

close all;
clear;

% basic parameters
d = 9;

to_save = 1;

[G, nodes_info] = hasse_diagram_dihedral(d);

num_subgroups = 0;
for j=1:size(G.Nodes,1)
    if nodes_info(j).is_dihedral
        num_subgroups = num_subgroups + 1;
    end
end
%num_subgroups = size(G.Nodes,1);
    
figure;
    
% Plot the graph without nodes
h = plot(G, 'Layout', 'layered', 'Direction', 'down', ...
    'NodeLabel', {}, 'Marker', 'none', 'LineWidth', 1.5);
highlight(h, 'Edges', 1:numedges(G), 'EdgeColor', 'k', 'LineWidth', 1.2);

% Add transparent circles for nodes
hold on;
scatter(h.XData, h.YData, 650, 'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', 'w', 'LineWidth', 1.5);

% Highlight dihedral nodes with bold circles
dihedral_idx = find([nodes_info.is_dihedral]);  % indices of dihedral nodes
dihedral_idx = and(mod(d,2*[(nodes_info.order)])==0,[nodes_info.is_dihedral]);
scatter(h.XData(dihedral_idx), h.YData(dihedral_idx), ...
    1000, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'none', 'LineWidth', 3);

hold off;

% Determine node colors: red for one entries in D, cyan otherwise (zero)
%node_colors = repmat([1 0 0], num_subgroups, 1);  % default cyan (RGB [0 1 1])
%zero_indices = find(D == 0);
%node_colors(zero_indices, :) = repmat([0 1 1], length(zero_indices), 1); % red
% h.NodeColor = node_colors;

% % Print subgroup labels inside nodes
% for i = 1:num_subgroups
%     x = h.XData(i);
%     y = h.YData(i);
%     text(x, y, sprintf('D_{%d,j}', nodes_info(i).order), ...
%         'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 10, 'FontWeight', 'bold');
% end
% for i = (num_subgroups+1):length(nodes_info)
%     x = h.XData(i);
%     y = h.YData(i);
%     text(x, y, sprintf('C_{%d}', nodes_info(i).order), ...
%         'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 10, 'FontWeight', 'bold');
% end


% Print subgroup labels inside nodes
for i = 1:length(nodes_info)
    x = h.XData(i);
    y = h.YData(i);

    if nodes_info(i).is_dihedral
        % Inside label (bold D)
        text(x, y, sprintf('D_{%d}', nodes_info(i).order), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
            'FontSize', 10, 'FontWeight', 'bold');

        % Outside label (j=0,...,d/m-1)
        m = nodes_info(i).order;
        j_max = d/m - 1;  % integer since m | d
        if j_max>0
            if j_max>1
                label_str = sprintf('j = 0,1,...,%d', j_max);
            else
                label_str = 'j = 0,1';
            end

        % Place slightly above the node
       % text(x+.6, y, label_str, ...
       %     'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
       %     'FontSize', 12, 'FontAngle', 'italic', 'Color', 'b');
        end
    else
        % Inside label for cyclic nodes
        text(x, y, sprintf('C_{%d}', nodes_info(i).order), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
            'FontSize', 10, 'FontWeight', 'bold');
    end
end

%place_dihedral_labels_lr(h, G, nodes_info, d, ...
 %   'BaseOffset', 0.35, ...      % base horizontal offset from node center
 %   'NodeClearance', 0.1, ...   % min dist from any node center
 %   'LabelClearance', 0.02, ...  % padding around label bbox
 %   'EdgeClearance', 0.02);      % inflate bbox to steer clear of edges
% After your scatter / labels inside nodes, call:

%place_dihedral_labels_lr(h, G, nodes_info, d, ...
%  'FontSize', 12, 'FontWeight','bold', 'Color','k', ...
%  'ShowBackground', true, 'BackgroundColor',[1 1 1], 'BackgroundMargin', 3, ...
%  'Halo', true, 'HaloColor',[1 1 1], 'HaloWidthPx', 3, ...
%  'Scales', [1 1.5 2 3 4 5]);

if to_save
    folder_name = ['DAG_dihedral_',num2str(d)];
    if isfolder(folder_name)
        cd(folder_name);
    else
        mkdir(folder_name)
        cd(folder_name)
    end

    name_it = ['hasse_dihedral_d_', num2str(d)] ;
    f = gcf;
    f.Position = [100 100 600 500];
    saveas(f, name_it ,'fig');
    saveas(f,name_it,'jpg');
    exportgraphics(f, [name_it,'.pdf'], 'ContentType', 'vector');
%   save('data_for_run');
    cd '../'
end
