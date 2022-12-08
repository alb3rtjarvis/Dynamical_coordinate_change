function [line] = line_between_points(pts,dx)
% Draws lines between pts with step size dt
if ~exist('dx','var')
    dx = 1e-5;
end
if size(pts,2)>2
    pts = pts';
end
line = cell(1,size(pts,1));
for k = 1:size(pts,1)-1
    if abs(pts(k,1)-pts(k+1,1))>1e-8
        if pts(k,1)<pts(k+1,1)
            p1 = pts(k,:); p2 = pts(k+1,:);
        else
            p1 = pts(k+1,:); p2 = pts(k,:);
        end
        m = slope(p1,p2);
        b = p1(2)-m*p1(1);
        x = p1(1):dx:p2(1);
        y = m*x+b; 
    else
        if pts(k,2)<pts(k+1,2)
            p1 = pts(k,:); p2 = pts(k+1,:);
        else
            p1 = pts(k+1,:); p2 = pts(k,:);
        end
        y = p1(2):dx:p2(2);
        x = p1(1)*ones(size(y));
    end
    line{k} = [x',y'];
end
if abs(pts(end,1)-pts(1,1))>1e-8
    if pts(end,1)<pts(1,1)
        p1 = pts(end,:); p2 = pts(1,:);
    else
        p1 = pts(1,:); p2 = pts(end,:);
    end
    m = slope(p1,p2);
    b = p1(2)-m*p1(1);
    x = p1(1):dx:p2(1);
    y = m*x+b; 
else
    if pts(end,2)<pts(1,2)
        p1 = pts(end,:); p2 = pts(1,:);
    else
        p1 = pts(1,:); p2 = pts(end,:);
    end
    y = p1(2):dx:p2(2);
    x = p1(1)*ones(size(y));
end
line{end} = [x',y'];
end
function m = slope(p1,p2)
    m = (p2(2)-p1(2))/(p2(1)-p1(1));
end

