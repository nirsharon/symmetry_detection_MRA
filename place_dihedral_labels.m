function warnings = place_dihedral_labels(h, G, nodes_info, d, varargin)
%PLACE_DIHEDRAL_LABELS  Add outside labels for dihedral nodes, avoiding nodes, labels, and edges.
%
% warnings = place_dihedral_labels(h, G, nodes_info, d, 'Name', value, ...)
%
% Required:
%   h           - GraphPlot object returned by plot(G, ...)
%   G           - graph object
%   nodes_info  - struct array with fields: is_dihedral (logical), order (integer)
%   d           - base group order used to compute j range: j = 0..(d/m - 1)
%
% Name-Value options:
%   'BaseOffset'     (0.12)  base offset (data units) from node center
%   'NodeClearance'  (0.10)  min distance (data units) from any node center
%   'LabelClearance' (0.02)  padding around label bbox (data units)
%   'EdgeClearance'  (0.02)  extra inflation of label bbox to avoid edges
%   'FontSize'       (9)     external label font size
%   'FontAngle'      ('italic')
%   'Color'          ('b')
%   'Candidates'     (8 dirs) candidate directions to try (rows of [dx dy])
%
% Returns:
%   warnings - cell array of warning strings (fallback placements, etc.)

% ---- Parse inputs ----
p = inputParser;
p.addParameter('BaseOffset', 0.22, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('NodeClearance', 0.10, @(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter('LabelClearance', 0.02, @(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter('EdgeClearance', 0.02, @(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter('FontSize', 12, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('FontAngle', 'italic', @(x)ischar(x)||isstring(x));
p.addParameter('Color', 'b');
p.addParameter('Candidates', [ 0  1; 1  0; 0 -1; -1  0; 1  1; -1  1; 1 -1; -1 -1 ], ...
    @(x)isnumeric(x)&&size(x,2)==2);
p.parse(varargin{:});
opt = p.Results;

ax = ancestor(h, 'axes');
hold(ax, 'on');

X = h.XData(:); Y = h.YData(:);
dihedral_idx = find([nodes_info.is_dihedral]);

% Precompute edges as segment endpoints in data coordinates
E = G.Edges.EndNodes;
sIdx = E(:,1);   % source node indices
tIdx = E(:,2);   % target node indices
% Fallback if nodes are numbered 1..n and Names not used:
if any(sIdx==0) || any(tIdx==0)
    try
        sIdx = str2double(E(:,1));
        tIdx = str2double(E(:,2));
    catch
        % if still failing, assume numeric node IDs already aligned with X/Y
        sIdx = (1:size(E,1))'; %#ok<NASGU>
        tIdx = sIdx;
    end
end
edgeSegs = [X(sIdx), Y(sIdx), X(tIdx), Y(tIdx)]; % [x1 y1 x2 y2] per edge

placed_extents = [];   % store visible label bboxes [x y w h]
warnings = {};

for k = 1:numel(dihedral_idx)
    i = dihedral_idx(k);
    x = X(i); y = Y(i);
    m = nodes_info(i).order;
    jmax = d/m - 1;
    label_str = sprintf('j = 0,1,...,%d', jmax);

    placed = false;
    for c = 1:size(opt.Candidates,1)
        dx = opt.Candidates(c,1) * opt.BaseOffset;
        dy = opt.Candidates(c,2) * opt.BaseOffset;

        % Measure bbox with invisible text at candidate position
        th = text(ax, x+dx, y+dy, label_str, ...
            'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
            'FontSize', opt.FontSize, 'FontAngle', opt.FontAngle, ...
            'Color', opt.Color, 'Units','data', 'Visible','off', 'Clipping','on');
        ext = get(th, 'Extent'); % [x y w h] in data units
        delete(th);

        % Pad bbox
        ext_pad = [ext(1)-opt.LabelClearance, ext(2)-opt.LabelClearance, ...
                   ext(3)+2*opt.LabelClearance, ext(4)+2*opt.LabelClearance];

        % 1) Check node clearance: bbox must be at least NodeClearance away from any node center
        if ~bbox_clears_nodes(ext_pad, X, Y, opt.NodeClearance)
            continue
        end

        % 2) Check label collision: bbox must not overlap previously placed bboxes
        if ~bbox_clears_labels(ext_pad, placed_extents)
            continue
        end

        % 3) Check edge clearance: bbox (inflated) must not intersect any edge segment
        if ~bbox_clears_edges(ext_pad, edgeSegs, opt.EdgeClearance)
            continue
        end

        % Place visible text
        th2 = text(ax, x+dx, y+dy, label_str, ...
            'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
            'FontSize', opt.FontSize, 'FontAngle', opt.FontAngle, ...
            'Color', opt.Color, 'Units','data', 'Clipping','on');
        placed_extents = [placed_extents; get(th2,'Extent')]; %#ok<AGROW>
        placed = true;
        break
    end

    if ~placed
        % Fallback: above the node (still respecting style; may overlap)
        th2 = text(ax, x, y+opt.BaseOffset, label_str, ...
            'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
            'FontSize', opt.FontSize, 'FontAngle', opt.FontAngle, ...
            'Color', opt.Color, 'Units','data', 'Clipping','on');
        placed_extents = [placed_extents; get(th2,'Extent')]; %#ok<AGROW>
        warnings{end+1} = sprintf('Label near node %d placed with fallback; may overlap.', i); %#ok<AGROW>
    end
end

% Optional console summary
if ~isempty(warnings)
    fprintf('[place_dihedral_labels] %d fallback(s):\n', numel(warnings));
    fprintf('  - %s\n', warnings{:});
end
end

% ---------- Helpers ----------
function ok = bbox_clears_nodes(ext, X, Y, minDist)
% Distance from point to rectangle; require >= minDist for all nodes
x1 = ext(1); y1 = ext(2); x2 = x1 + ext(3); y2 = y1 + ext(4);
px = max(x1, min(X, x2));
py = max(y1, min(Y, y2));
dist = hypot(X - px, Y - py);
ok = all(dist >= minDist);
end

function ok = bbox_clears_labels(ext, placed_extents)
% No overlap with any previously placed bbox
if isempty(placed_extents)
    ok = true; return
end
r1 = ext;
noOverlap = (r1(1)+r1(3) < placed_extents(:,1)) | ...
            (placed_extents(:,1)+placed_extents(:,3) < r1(1)) | ...
            (r1(2)+r1(4) < placed_extents(:,2)) | ...
            (placed_extents(:,2)+placed_extents(:,4) < r1(2));
ok = all(noOverlap);
end

function ok = bbox_clears_edges(ext, edgeSegs, inflate)
% Reject if any edge segment intersects the (inflated) rectangle
r = [ext(1)-inflate, ext(2)-inflate, ext(3)+2*inflate, ext(4)+2*inflate];
for e = 1:size(edgeSegs,1)
    p1 = edgeSegs(e,1:2);
    p2 = edgeSegs(e,3:4);
    if seg_intersects_rect(p1, p2, r)
        ok = false; return
    end
end
ok = true;
end

function tf = seg_intersects_rect(p1, p2, rect)
% rect = [x y w h]
x1 = rect(1); y1 = rect(2); x2 = x1 + rect(3); y2 = y1 + rect(4);
% 1) Endpoint inside rect?
if point_in_rect(p1, rect) || point_in_rect(p2, rect)
    tf = true; return
end
% 2) Segment intersects any of the 4 rectangle edges?
R = [x1 y1; x2 y1; x2 y2; x1 y2];
edges = [R; R(1,:)]; % close loop
for k = 1:4
    q1 = edges(k,:); q2 = edges(k+1,:);
    if segments_intersect(p1, p2, q1, q2)
        tf = true; return
    end
end
tf = false;
end

function tf = point_in_rect(p, rect)
x = p(1); y = p(2);
x1 = rect(1); y1 = rect(2); x2 = x1 + rect(3); y2 = y1 + rect(4);
tf = (x >= x1) && (x <= x2) && (y >= y1) && (y <= y2);
end

function tf = segments_intersect(p1, p2, q1, q2)
% Standard orientation test for 2D segment intersection
o1 = orient(p1, p2, q1);
o2 = orient(p1, p2, q2);
o3 = orient(q1, q2, p1);
o4 = orient(q1, q2, p2);
if o1 ~= o2 && o3 ~= o4
    tf = true; return
end
% Colinear special cases
tf = on_seg(p1,q1,p2) || on_seg(p1,q2,p2) || on_seg(q1,p1,q2) || on_seg(q1,p2,q2);
end

function val = orient(a, b, c)
val = sign( (b(2)-a(2))*(c(1)-b(1)) - (b(1)-a(1))*(c(2)-b(2)) );
end

function tf = on_seg(a, b, c)
tf = min(a(1),c(1)) <= b(1) && b(1) <= max(a(1),c(1)) && ...
     min(a(2),c(2)) <= b(2) && b(2) <= max(a(2),c(2)) && ...
     orient(a,c,b)==0;
end
