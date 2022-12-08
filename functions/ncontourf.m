function [ch] = ncontourf(X,Y,Z,varargin)
%%% Shortcut contour plot that removes edgecolor and sets equal aspect
%%% ratio
[~,ch] = contourf(X,Y,Z,varargin{:});
set(ch,'edgecolor','none'); % remove the lines between contour levels
daspect([1 1 1]); % set the aspect ratio to 1:1
drawnow
end

