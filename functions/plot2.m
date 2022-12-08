function [] = plot2(x,varargin)
%%% Simplified function for plotting when all data is stored in one array
    plot(x(:,1),x(:,2),varargin{:})
end

