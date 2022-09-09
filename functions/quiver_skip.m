function [] = quiver_skip(X,Y,U,V,s,varargin)
%%% Function to produce quiver plot with skip of s size in data where
%%% vector field is defined by meshgrids U and V
xmin = min(min(X));
xmax = max(max(X));
ymin = min(min(Y));
ymax = max(max(Y));
if size(X,1) == 1 || size(X,2) == 1
    if size(U,1) == 1 || size(U,2) == 1
        quiver(X(1:s:end),Y(1:s:end),U(1:s:end),V(1:s:end),varargin{:});
    else
        quiver(X(1:s:end),Y(1:s:end),U(1:s:end,1:s:end),V(1:s:end,1:s:end),varargin{:});
    end
else
    if size(U,1) == 1 || size(U,2) == 1
        quiver(X(1:s:end,1:s:end),Y(1:s:end,1:s:end),U(1:s:end),V(1:s:end),varargin{:});
    else
        quiver(X(1:s:end,1:s:end),Y(1:s:end,1:s:end),U(1:s:end,1:s:end),V(1:s:end,1:s:end),varargin{:});
    end
end
xlim([xmin xmax]);
ylim([ymin ymax]);
drawnow
end

