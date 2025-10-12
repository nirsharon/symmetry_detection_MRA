function plot_hasse_clean(G, nodes_info, varargin)
% PLOT_HASSE_CLEAN  Deterministic Hasse layout with levels by subgroup order.
% Options:
%   'RotationsOnly', true/false   % if true, hides reflections & D_{2m} copies
%   'ShowLabels', true/false
%   'MinGap', 0.8                 % horizontal spacing scale

p = inputParser;
addParameter(p,'RotationsOnly',false,@islogical);
addParameter(p,'ShowLabels',true,@islogical);
addParameter(p,'MinGap',0.8,@(x)isnumeric(x)&&isscalar(x));
parse(p,varargin{:});
rotOnly   = p.Results.RotationsOnly;
showLbl   = p.Results.ShowLabels;
mingap    = p.Results.MinGap;

% ---- Optional filtering to mimic the cyclic lattice ----
keep = true(numel(nodes_info),1);
if rotOnly
    % Keep only 1, C_m, and the top D_{2n}; drop reflections & D_{2m} (m<n)
    for i=1:numel(nodes_info)
        t = getfield(nodes_info(i),'type'); %#ok<GFLD>
        if ~(strcmp(t,'1') || strcmp(t,'C') || strcmp(t,'G'))
            keep(i) = false;
        end
    end
    % induce subgraph
    G = subgraph(G, find(keep));
    nodes_info = nodes_info(keep);
end

% ---- Compute levels: larger-order higher ----
ord = arrayfun(@(s)s.order, nodes_info);
% Unique orders (descending) define levels
[uOrders,~,idxU] = unique(ord);
[~,ordPos] = sort(uOrders,'descend');
level = ordPos(idxU);   % 1 = top level

% Spread nodes horizontally within the same level by a stable key:
key = arrayfun(@(s) level_key(s), nodes_info, 'UniformOutput', false);
x = zeros(numel(nodes_info),1); y = zeros(numel(nodes_info),1);
for L = unique(level(:))'
    I = find(level==L);
    [~,perm] = sort(key(I));            % stable order inside the level
    k = numel(I);
    xs = ((1:k) - (k+1)/2) * mingap;    % centered positions
    x(I(perm)) = xs;
    y(I(perm)) = max(level)-L;          % top level gets largest y
end

% Prepare labels
if isfield(nodes_info,'label')
    labels = {nodes_info.label};
else
    % fallback
    labels = arrayfun(@(s)sprintf('(%d)',s.order), nodes_info, 'UniformOutput', false);
end
if ~showLbl
    labels = repmat({''}, size(labels));
end

% Plot with fixed positions
h = plot(G,'XData',x,'YData',y,'MarkerSize',7,'LineWidth',0.75);
set(h,'NodeLabel',labels);
set(gca,'Visible','off');

% Style by type (if available)
if isfield(nodes_info,'type')
    C = zeros(numel(nodes_info),3);
    for i=1:numel(nodes_info)
        switch nodes_info(i).type
            case '1', C(i,:) = [0 0 0];
            case 'G', C(i,:) = [0 0 0];
            case 'C', C(i,:) = [0.10 0.45 0.80];
            case {'D'}, C(i,:) = [0.80 0.40 0.10];
            case {'R','Rall'}, C(i,:) = [0.20 0.60 0.20];
            otherwise, C(i,:) = [0.5 0.5 0.5];
        end
    end
    h.NodeColor = C;
end

% A simple title
if any(strcmp({nodes_info.type},'G'))
    top = nodes_info(strcmp({nodes_info.type},'G'));
    if ~isempty(top)
        title(sprintf('Hasse diagram of %s', top(1).label), 'Interpreter','tex');
    end
end

% ------------ helpers ------------
function s = level_key(node)
    % Order first; then a stable secondary key by type and parameters, so
    % nodes within a level have predictable left-to-right placement.
    t = '(?)'; a=''; 
    if isfield(node,'type'), t=node.type; end
    if isfield(node,'params')
        P = node.params;
        if isfield(P,'m'), a = [a, sprintf('m=%d;',P.m)]; end %#ok<AGROW>
        if isfield(P,'a'), a = [a, sprintf('a=%d;',P.a)]; end %#ok<AGROW>
        if isfield(P,'k'), a = [a, sprintf('k=%d;',P.k)]; end %#ok<AGROW>
        if isfield(P,'copies'), a = [a, sprintf('c=%d;',P.copies)]; end %#ok<AGROW>
    end
    % type order: 1 < C < R < D < G inside the same order (rare but deterministic)
    tRank = struct('one',0,'C',1,'R',2,'Rall',2,'D',3,'G',4);
    if isfield(tRank, t), r = tRank.(t); else, r = 5; end
    s = sprintf('%05d|%02d|%s', -node.order, r, a); % negative so higher order sorts first
end
end
