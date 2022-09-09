function [] = qquiver(x,y,v)
dims = numel(size(v));
    if dims == 2
        quiver(x,y,v(:,1),v(:,2));
    else
        quiver(x,y,v(:,:,1),v(:,:,2));
    end
end
