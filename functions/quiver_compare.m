function [] = quiver_compare(x,y,v1,v2,varargin)
%%% Function to compare two quiver plots defined on x,y where the vector
%%% fields are given by v1 and v2, can be in meshgrid format or arrays
dims = numel(size(v1));
    if dims == 2
        quiver(x,y,v1(:,1),v1(:,2),varargin{:})
        hold on;
        quiver(x,y,v2(:,1),v2(:,2),varargin{:})
    else
        quiver(x,y,v1(:,:,1),v1(:,:,2),varargin{:})
        hold on;
        quiver(x,y,v2(:,:,1),v2(:,:,2),varargin{:})
    end
end
