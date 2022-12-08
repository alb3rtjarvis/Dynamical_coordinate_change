function [] = scatter2(xyvec,varargin)
%%% Shortcut scatter plot for when data is all stored in xyvec
if size(xyvec,2)==2
    scatter(xyvec(:,1),xyvec(:,2),varargin{:});
else
    scatter(xyvec(1,:),xyvec(2,:),varargin{:});
end

