function warnings = place_dihedral_labels_lr(h, G, nodes_info, d, varargin)
% Outside labels for dihedral nodes, placed LEFT/RIGHT with adaptive offset.
% Avoids nodes, already-placed labels, and edges. Assumes numeric node indices.

% ---- Options ----
p = inputParser;
p.addParameter('BaseOffset', 0.33);
p.addParameter('NodeClearance', 0.10);
p.addParameter('LabelClearance', 0.02);
p.addParameter('EdgeClearance', 0.02);
p.addParameter('FontSize', 10);                 % ⬅ bigger default
p.addParameter('FontAngle', 'italic');
p.addParameter('FontWeight', 'bold');           % ⬅ bold by default
p.addParameter('Color', 'k');                   % ⬅ high contrast
p.addParameter('Scales', [1 1.5 2 3]);        % ⬅ search further if tight
% Visibility helpers:
p.addParameter('ShowBackground', true);         % ⬅ white bg behind text
p.addParameter('BackgroundColor', [1 1 1]);     % solid white (no alpha)
p.addParameter('BackgroundMargin', 2);          % pixels
p.addParameter('EdgeColor', 'none');            % bg border
p.addParameter('Halo', true);                   % ⬅ draw a subtle halo
p.addParameter('HaloColor', [1 1 1]);           % white halo
p.addParameter('HaloWidthPx', 2);               % halo thickness (px)
p.parse(varargin{:});
o = p.Results;

ax = ancestor(h,'axes'); hold(ax,'on');
X = h.XData(:); Y = h.YData(:);

% Dihedral nodes
dihedral_idx = find([nodes_info.is_dihedral]);

% Edges -> line segments (numeric)
E = G.Edges.EndNodes;
sIdx = E(:,1); tIdx = E(:,2);
edgeSegs = [X(sIdx), Y(sIdx), X(tIdx), Y(tIdx)];

% Track placed label bounding boxes
placed_extents = [];
warnings = {};

for k = 1:numel(dihedral_idx)
    i = dihedral_idx(k);
    x = X(i); y = Y(i);
    m = nodes_info(i).order;
    jmax = d/m - 1;

    % ---- Your label rule ----
    if jmax > 2
        label_str = sprintf('j = 0,1,...,%d', jmax);
    else
        switch jmax
            case 2
                label_str = 'j = 0,1,2';
            case 1
                label_str = 'j = 0,1';
            otherwise
                label_str = '';
        end
    end
    if isempty(label_str)
        continue; % skip if no label
    end

    % Prefer side with more axis room
    xlim = ax.XLim;
    roomRight = xlim(2) - x;
    roomLeft  = x - xlim(1);
    if roomLeft > roomRight
        candidates = [-1 0; 1 0]; % left first
    else
        candidates = [ 1 0; -1 0]; % right first
    end

    placed = false;
    for s = 1:numel(o.Scales)
        scale = o.Scales(s);
        for c = 1:size(candidates,1)
            dx = candidates(c,1) * o.BaseOffset * scale;
            dy = 0;

            % Invisible text to measure bbox (no bg here for tight extent)
            th = text(ax, x+dx, y+dy, label_str, ...
                'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
                'FontSize', o.FontSize, 'FontAngle', o.FontAngle, ...
                'FontWeight', o.FontWeight, 'Units','data', 'Visible','off');
            ext = get(th,'Extent'); delete(th);

            % Pad bbox
            ext_pad = [ext(1)-o.LabelClearance, ext(2)-o.LabelClearance, ...
                       ext(3)+2*o.LabelClearance, ext(4)+2*o.LabelClearance];

            % Checks: nodes, labels, edges
            if ~bbox_clears_nodes(ext_pad, X, Y, o.NodeClearance), continue; end
            if ~bbox_clears_labels(ext_pad, placed_extents), continue; end
            if ~bbox_clears_edges(ext_pad, edgeSegs, o.EdgeClearance), continue; end

            % ---- Place visible label (with visibility helpers) ----
            th2 = place_text_visible(ax, x+dx, y+dy, label_str, o);
            uistack(th2,'top');                               % ensure on top
            placed_extents = [placed_extents; get(th2,'Extent')]; %#ok<AGROW>
            placed = true;
            break
        end
        if placed, break; end
    end

    if ~placed
        % Fallback: base offset to preferred side
        dx = candidates(1,1) * o.BaseOffset;
        th2 = place_text_visible(ax, x+dx, y, label_str, o);
        uistack(th2,'top');
        placed_extents = [placed_extents; get(th2,'Extent')]; %#ok<AGROW>
        warnings{end+1} = sprintf('Label near node %d placed with fallback.', i); %#ok<AGROW>
    end
end

if ~isempty(warnings)
    fprintf('[place_dihedral_labels_lr] %d fallback(s):\n', numel(warnings));
    fprintf('  - %s\n', warnings{:});
end
end

% ---- Helpers ----
function th = place_text_visible(ax, x, y, str, o)
% Draws a halo (optional), then the text with background (optional), returns text handle.

% Convert halo width (px) to data units:
dx_data = px2data_x(ax, o.HaloWidthPx);
dy_data = px2data_y(ax, o.HaloWidthPx);

% Halo: draw 8 surrounding copies in halo color behind the main text
if o.Halo
    dirs = [ 1 0; -1 0; 0 1; 0 -1;  1 1; 1 -1; -1 1; -1 -1 ];
    for ii = 1:size(dirs,1)
        text(ax, x + dirs(ii,1)*dx_data, y + dirs(ii,2)*dy_data, str, ...
            'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
            'FontSize', o.FontSize, 'FontAngle', o.FontAngle, ...
            'FontWeight', o.FontWeight, 'Color', o.HaloColor, ...
            'Units','data', 'Clipping','on', 'HitTest','off');
    end
end

% Main text (with optional background box)
th = text(ax, x, y, str, ...
    'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
    'FontSize', o.FontSize, 'FontAngle', o.FontAngle, ...
    'FontWeight', o.FontWeight, 'Color', o.Color, ...
    'Units','data', 'Clipping','on', 'HitTest','off');

if o.ShowBackground
    set(th, 'BackgroundColor', o.BackgroundColor, ...
            'Margin', o.BackgroundMargin, ...
            'EdgeColor', o.EdgeColor);
end
end

function v = px2data_x(ax, px)
% Convert pixels to data units along x
oldUnits = ax.Units; ax.Units = 'pixels';
pos = ax.Position; ax.Units = oldUnits;
v = (ax.XLim(2)-ax.XLim(1)) * (px / max(pos(3),1));
end

function v = px2data_y(ax, px)
% Convert pixels to data units along y
oldUnits = ax.Units; ax.Units = 'pixels';
pos = ax.Position; ax.Units = oldUnits;
v = (ax.YLim(2)-ax.YLim(1)) * (px / max(pos(4),1));
end

function ok = bbox_clears_nodes(ext, X, Y, minDist)
x1=ext(1); y1=ext(2); x2=x1+ext(3); y2=y1+ext(4);
px = max(x1, min(X, x2));
py = max(y1, min(Y, y2));
dist = hypot(X - px, Y - py);
ok = all(dist >= minDist);
end

function ok = bbox_clears_labels(ext, placed_extents)
if isempty(placed_extents), ok = true; return; end
r1 = ext;
noOverlap = (r1(1)+r1(3) < placed_extents(:,1)) | ...
            (placed_extents(:,1)+placed_extents(:,3) < r1(1)) | ...
            (r1(2)+r1(4) < placed_extents(:,2)) | ...
            (placed_extents(:,2)+placed_extents(:,4) < r1(2));
ok = all(noOverlap);
end

function ok = bbox_clears_edges(ext, edgeSegs, inflate)
r = [ext(1)-inflate, ext(2)-inflate, ext(3)+2*inflate, ext(4)+2*inflate];
for e = 1:size(edgeSegs,1)
    p1 = edgeSegs(e,1:2); p2 = edgeSegs(e,3:4);
    if seg_intersects_rect(p1, p2, r), ok = false; return; end
end
ok = true;
end

function tf = seg_intersects_rect(p1,p2,rect)
x1=rect(1); y1=rect(2); x2=x1+rect(3); y2=y1+rect(4);
if point_in_rect(p1,rect) || point_in_rect(p2,rect), tf = true; return; end
R=[x1 y1; x2 y1; x2 y2; x1 y2]; edges=[R;R(1,:)];
for k=1:4
    if segments_intersect(p1,p2,edges(k,:),edges(k+1,:)), tf=true; return; end
end
tf=false;
end

function tf = point_in_rect(p,rect)
x=p(1); y=p(2); x1=rect(1); y1=rect(2); x2=x1+rect(3); y2=y1+rect(4);
tf=(x>=x1)&&(x<=x2)&&(y>=y1)&&(y<=y2);
end

function tf = segments_intersect(p1,p2,q1,q2)
o1=orient(p1,p2,q1); o2=orient(p1,p2,q2);
o3=orient(q1,q2,p1); o4=orient(q1,q2,p2);
if o1~=o2 && o3~=o4, tf=true; return; end
tf=on_seg(p1,q1,p2)||on_seg(p1,q2,p2)||on_seg(q1,p1,q2)||on_seg(q1,p2,q2);
end

function val = orient(a,b,c)
val = sign((b(2)-a(2))*(c(1)-b(1))-(b(1)-a(1))*(c(2)-b(2)));
end

function tf = on_seg(a,b,c)
tf = min(a(1),c(1))<=b(1)&&b(1)<=max(a(1),c(1)) && ...
     min(a(2),c(2))<=b(2)&&b(2)<=max(a(2),c(2)) && ...
     orient(a,c,b)==0;
end
