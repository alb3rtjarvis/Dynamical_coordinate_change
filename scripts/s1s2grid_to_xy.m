% Starting from a predefined grid in s1s2, xy coordinates are computed by
% integrating tensorlines of the right Cauchy-Green strain tensor
%% Define xygrid, parameters, vfield and eigenvector fields on xy points for plotting purposes
clear; clc;
addpath('../functions')
set(0,'defaultAxesFontSize',24);
set(0, 'DefaultLineLineWidth', 2);
set(0,'defaulttextInterpreter','latex');
nx = 101;
ny = 51;
x = linspace(0,2,nx);
y = linspace(0,1,ny);
[X,Y] = meshgrid(x,y);
t0 = 0;
T = 10;
tspan = linspace(t0,T,501);
vfield = @double_gyre;
[l_max,l_min,ev_max,ev_min] = CG_ev_field(vfield,tspan,x,y,T);
%% Create functions for eigenvector fields to be used in integration

% Define domain
figure
qquiver(X,Y,ev_min)
hold on


dX = 0.001;
x_cut = [0.23,1.95];
x_left = x_cut(1)*ones(size(y));
x_right = x_cut(2)*ones(size(y));
y_cut = [0.05,0.95];
y_bot = y_cut(1)*ones(size(x));
y_top = y_cut(2)*ones(size(x));
plot(x_left,y,x_right,y)
plot(x,y_bot,x,y_top)

d1x = 0.661; d1y = [0.539,0.641];
d2x = 1.439; d2y = [0.359,0.501];
d1x_line = 0:dX:d1x;
d1y_line_bot = d1y(1)*ones(size(d1x_line));
d1y_line_top = d1y(2)*ones(size(d1x_line));

d2x_line = d2x:dX:x(end);
d2y_line_bot = d2y(1)*ones(size(d2x_line));
d2y_line_top = d2y(2)*ones(size(d2x_line));

plot(d1x_line,d1y_line_top,d1x_line,d1y_line_bot);
plot(d2x_line,d2y_line_top,d2x_line,d2y_line_bot);

d1vy = d1y(1):dX:d1y(2);
d1vx = d1x*ones(size(d1vy));
d2vy = d2y(1):dX:d2y(2);
d2vx = d2x*ones(size(d2vy));
plot(d1vx,d1vy,d2vx,d2vy)
%% Create struct containing boundary lines, for bound_struct in integration
bs = struct();
bs(1).x = d1x_line';
bs(1).y = d1y_line_top';
bs(2).x = d1x_line';
bs(2).y = d1y_line_bot';
bs(3).x = d2x_line';
bs(3).y = d2y_line_top';
bs(4).x = d2x_line';
bs(4).y = d2y_line_bot';
bs(5).x = d1vx';
bs(5).y = d1vy';
bs(6).x = d2vx';
bs(6).y = d2vy';
bs(7).x = x_left';
bs(7).y = y';
bs(8).x = x';
bs(8).y = y_top';
%% Create functions for eigenvector fields to be used in integration


% Define domain and origin for arc length axis
ub = [max(x),max(y)];
lb = [min(x),min(y)];
D = [lb;ub];
O = [1.098,0.5242];

evfield1 = @(t,x) CG_eigenvector(t0,x,@double_gyre,T);
evfield2 = @(t,x) CG_eigenvector(t0,x,@double_gyre,T,1);

% Integrate initial arclength axis tensorlines to determine maxlength
% before reaching degenerate points, also used for stop_line
[X_out1f,sfspan] = integrate_tensorlines_grow_bound(evfield1,0:.01:0.5,t0,O,D,T,bs);
[X_out1b,sbspan] = integrate_tensorlines_grow_bound(evfield1,0:-.01:-0.5,t0,O,D,T,bs);
hold on
plot2(X_out1f)
plot2(X_out1b)
%% Determine max length of arclength axis before reaching degenerate points
a1f = abs(sfspan(end));
a1b = abs(sbspan(end));
len_a1f = (floor(a1f*100))/100;
len_a1b = (floor(a1b*100))/100;
len = min(len_a1f,len_a1b);
% step size in arclength
da1 = 1e-2;
pts1 = 0:da1:len;
pts = flip(unique([flip(pts1),-pts1]));
s1axis_t = unique([flipud(X_out1f);X_out1b],'rows');
s1pts = flip(unique([flip(sfspan),sbspan]));
% Create interpolant for Stop_line in integration
s1axis = interp1(s1pts,s1axis_t,pts);
%% Perfrom coordinate transformation from arclength grid
fwd_span = 0:1e-3:2;
back_span = -fwd_span;

for k = 1:numel(pts)
    [s2fwd{k},s2f{k}] = integrate_tensorlines_grow_bound(evfield2,fwd_span,t0,s1axis(k,:),D,T,bs);
    [s2back{k},s2b{k}] = integrate_tensorlines_grow_bound(evfield2,back_span,t0,s1axis(k,:),D,T,bs);
end
%% Plot xy points corresponding to s1s2 grid
figure
qquiver(X,Y,ev_min)
hold on

plot(d1x_line,d1y_line_top,d1x_line,d1y_line_bot);
plot(d2x_line,d2y_line_top,d2x_line,d2y_line_bot);
plot(d1vx,d1vy,d2vx,d2vy)
plot2(s1axis)
for k = 1:numel(pts)
    plot2(s2fwd{k})
    plot2(s2back{k})
    drawnow
end

%% Clean data to make plotting easier
clear s2fc s2bc s2fwd_c s2back_c
n = 471;
i = 1;
idx = 1:numel(s2f);
for k = 1:numel(s2f)
    if numel(s2f{k}) < n || numel(s2b{k}) < n
        didx(i) = k;
        i = i+1;
    end
end
kidx = setdiff(idx,didx);
kn = numel(kidx);
kpts = pts(kidx);
ks1axis = s1axis(kidx,:);
i = 1;
for k = kidx
    s2fwd_c{i} = s2fwd{k}(1:10:n,:);
    s2fc{i} = s2f{k}(1:10:n)';
    s2back_c{i} = s2back{k}(1:10:n,:);
    s2bc{i} = s2b{k}(1:10:n)';
    i = i+1;
end
an = numel(s2bc{1});
%% Define arrays containing corresponding xy and s1s2 points
xyvec = [];
svec = [];
for k = 1:kn
    xyvec = [xyvec;flipud(s2back_c{k});s2fwd_c{k}(2:end,:)];
    svec = [svec;ones(an,1)*kpts(k),flipud(s2bc{k});ones(an-1,1)*kpts(k),s2fc{k}(2:end,:)];    
end
%%

%%
% %%
% 
% %%
figure
qquiver(X,Y,ev_min)
hold on
% plot(d1x_line,d1y_line_top,d1x_line,d1y_line_bot);
% plot(d2x_line,d2y_line_top,d2x_line,d2y_line_bot);
% plot(d1vx,d1vy,d2vx,d2vy)
plot2(ks1axis)
for k = 1:numel(s2fwd_c)
    plot2(s2fwd_c{k})
    plot2(s2back_c{k})
    drawnow
end
    
