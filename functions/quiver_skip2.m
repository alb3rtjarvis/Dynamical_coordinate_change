function [] = quiver_skip2(X,Y,VF,s,varargin)
xmin = min(min(X));
xmax = max(max(X));
ymin = min(min(Y));
ymax = max(max(Y));
vfdim = numel(size(VF));
if size(X,1) == 1 || size(X,2) == 1
    if vfdim == 2
        quiver(X(1:s:end),Y(1:s:end),VF(1:s:end,1),VF(1:s:end,2),varargin{:});
    else
        quiver(X(1:s:end),Y(1:s:end),VF(1:s:end,1:s:end,1),VF(1:s:end,1:s:end,2),varargin{:});
    end
else
    if vfdim == 2
        quiver(X(1:s:end,1:s:end),Y(1:s:end,1:s:end),VF(1:s:end,1),VF(1:s:end,2),varargin{:});
    else
        quiver(X(1:s:end,1:s:end),Y(1:s:end,1:s:end),VF(1:s:end,1:s:end,1),VF(1:s:end,1:s:end,2),varargin{:});
    end
end
xlim([xmin xmax]);
ylim([ymin ymax]);
drawnow
end

