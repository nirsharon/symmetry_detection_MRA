function [] = draw_points_unit_circle(args)
%

pts = zeros(numel(args),2);
for j=1:numel(args)
    pts(j,1) = cos(args(j)); 
    pts(j,2) = sin(args(j)); 
end

figure;
scatter(pts(:,1),pts(:,2),'filled');