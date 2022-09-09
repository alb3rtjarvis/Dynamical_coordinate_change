% Starting from a predefined grid in xy, s1s2 coordinates are computed by
% integrating tensorlines of the right Cauchy-Green strain tensor
%% Define xygrid, parameters, vfield
clear; clc;
addpath('../functions')
set(0,'defaultAxesFontSize',24);
set(0, 'DefaultLineLineWidth', 2);
set(0,'defaulttextInterpreter','latex');
verify_plots = false;
nx = 101;
ny = 51;
x = linspace(0,2,nx);
y = linspace(0,1,ny);
[X,Y] = meshgrid(x,y);
t0 = 0;
T = 10;
intspan = [t0,T];
tspan = linspace(t0,T,501);
vfield = @double_gyre;
%% Create functions for eigenvector fields to be used in integration

% Define domain and origin for arc length axis
ub = [max(x),max(y)];
lb = [min(x),min(y)];
O = [1.098,0.5242];
D = [lb;ub];
evfield1 = @(t,x) CG_eigenvector(t0,x,vfield,T);
evfield2 = @(t,x) CG_eigenvector(t0,x,vfield,T,1);
%% Identify degenerate points in tensor field and define cut domains
x1 = [0.5,0.5];
x2 = [1.5,0.5];
options = optimset('TolFun',1e-13);
optfun = @(x) CG_search(t0,x,vfield,T);
[x1d] = fminsearch(optfun,x1,options);
[x2d] = fminsearch(optfun,x2,options);
X1d_aux = aux_grid8(x1d,5e-2);
X2d_aux = aux_grid8(x2d,5e-2);
X1d = zeros(size(X1d_aux));
X1d(1,:) = [];
X2d = size(X1d);
for k = 1:size(X1d,1)
    X1d(k,:) = fminsearch(optfun,X1d_aux(k+1,:),options);
    X2d(k,:) = fminsearch(optfun,X2d_aux(k+1,:),options);
end

% Identify the unique degenerate pairs and define box centers
degen_pair1 = uniquetol(X1d,1e-5,'ByRows',true);
degen_pair2 = uniquetol(X2d,1e-5,'ByRows',true);
degen_array = [degen_pair1;degen_pair2];
%% Integrate to obtain arclength axis   
X_out1f = integrate_tensorlines_grow(evfield1,0:.005:0.65,O,D,degen_array);
X_out1b = integrate_tensorlines_grow(evfield1,0:-.005:-0.65,O,D,degen_array);
%% Create circles enclosing degenerate points
dc = 5e-3;
circ_center1 = sum(degen_pair1,1)/2;
circ_center2 = sum(degen_pair2,1)/2;
circ_size = 0.8;
delta1 = circ_size*norm(degen_pair1(1,:)-degen_pair1(2,:));
delta2 = circ_size*norm(degen_pair2(1,:)-degen_pair2(2,:));
circ1 = flipud(aux_circle(circ_center1,delta1,dc)');
circ2 = flipud(aux_circle(circ_center2,delta2,dc)');
%% Create polygon areas which seperate different domains
% Create outer boundary for flag conditions in integrator
edge_tol = 9e-2;
left_edge = min(x)+edge_tol; right_edge = max(x)-edge_tol;
bot_edge = min(y)+edge_tol; top_edge = max(y)-edge_tol;
dp = dc*min(delta1,delta2);

% Create cuts to remove degenerate domains
cut_width = 1e-2;
cut1y = [circ_center1(2)+cut_width/2;circ_center1(2)-cut_width/2];
cut2y = [circ_center2(2)+cut_width/2;circ_center2(2)-cut_width/2];

% Find where cuts will intersect circles
[~,iy1t] = mink(abs(circ1(:,2)-cut1y(1)),2);
xcut1t = min(circ1(iy1t,:));

[~,iy1b] = mink(abs(circ1(:,2)-cut1y(2)),2);
xcut1b = min(circ1(iy1b,:));

[~,iy2t] = mink(abs(circ2(:,2)-cut2y(1)),2);
xcut2t = max(circ2(iy2t,:));

[~,iy2b] = mink(abs(circ2(:,2)-cut2y(1)),2);
xcut2b = max(circ2(iy2b,:));

% Create horizontal cuts
cut1xt = min(x):dp:xcut1t;
cut1yt = cut1y(1)*ones(size(cut1xt));
cut1t = [cut1xt;cut1yt];
cut1xb = min(x):dp:xcut1b;
cut1yb = cut1y(2)*ones(size(cut1xb));
cut1b = [cut1xb;cut1yb];

cut2xt = xcut2t:dp:max(x);
cut2yt = cut2y(1)*ones(size(cut2xt));
cut2t = [cut2xt;cut2yt];
cut2xb = xcut2b:dp:max(x);
cut2yb = cut2y(2)*ones(size(cut2xb));
cut2b = [cut2xb;cut2yb];

% Seperate circles into top and bottom pieces for integration direction
circ1t = circ1(circ1(:,2)>circ_center1(2),:);
circ1b = circ1(circ1(:,2)<circ_center1(2),:);
[~,order1t] = sort(circ1t(:,1));
circ1t = circ1t(order1t,:);
[~,order1b] = sort(circ1b(:,1));
circ1b = circ1b(order1b,:);

circ2t = circ2(circ2(:,2)>circ_center2(2),:);
circ2b = circ2(circ2(:,2)<circ_center2(2),:);
[~,order2t] = sort(circ2t(:,1));
circ2t = circ2t(order2t,:);
[~,order2b] = sort(circ2b(:,1));
circ2b = circ2b(order2b,:);

circ_skip1 = ceil(cut_width/(1.5*dc*delta1));
circ_skip2 = ceil(cut_width/(1.5*dc*delta2));
cut_skip1 = ceil(0.005*length(cut1t));
cut_skip2 = ceil(0.005*length(cut2t));
% Create top, bottom, and full branches
top_branch1 = [cut1t(:,1:end-cut_skip1-1)';circ1t(circ_skip1:end,:)];
bot_branch1 = [cut1b(:,1:end-cut_skip1-1)';circ1b(circ_skip1:end,:)];
full_branch1 = [top_branch1;flipud(bot_branch1)];
top_branch2 = [circ2t(1:end-circ_skip2-1,:);cut2t(:,cut_skip2:end)'];
bot_branch2 = [circ2b(1:end-circ_skip2-1,:);cut2b(:,cut_skip2:end)'];
full_branch2 = [top_branch2;flipud(bot_branch2)];
%% Create interpolant of arclength axis, for Stop_line in integration
ni = 100001;
s1_line = [flipud(X_out1f);X_out1b];
s1_line = unique(s1_line,'rows');
t = linspace(0,1,ni); 
%% Compute arc length for each point in axis
s1_arc = interparc(t,s1_line(:,1),s1_line(:,2));
s1_O = abs(s1_arc-O);
[r,~] = find(vecnorm(s1_O,2,2) == min(vecnorm(s1_O,2,2)));
for k = 1:ni
    if k < r
        s1_arc(k,3) = arclength(s1_arc(k:r,1),s1_arc(k:r,2));
    elseif k > r
        s1_arc(k,3) = -arclength(s1_arc(r:k,1),s1_arc(r:k,2));
    else
        s1_arc(k,3) = 0;
    end
end

%% Create polygons which will define xy points to be tranformed

polygon(1).x = [left_edge;left_edge;right_edge;right_edge];
polygon(1).y = [bot_edge;top_edge;top_edge;bot_edge];
polygon(2).x = full_branch1(:,1);
polygon(2).y = full_branch1(:,2);
polygon(3).x = full_branch2(:,1);
polygon(3).y = full_branch2(:,2);
xq = X(:); 
yq = Y(:);
in = false(numel(xq),3);
% Identify which xy points lie within domain to be transformed
for k = 1:numel(polygon)
    in(:,k) = inpolygon(xq,yq,polygon(k).x,polygon(k).y);
end

for k = 1:numel(polygon)-1
    in(in(:,k+1),1) = 0;
end
xt = xq(in(:,1)); yt = yq(in(:,1));  
X0 = [xt,yt];
%% Define regions of the domain corresponding to init integration direction

n2 = 10001;
t2 = linspace(0,1,n2);
s1_line2 = interparc(t2,s1_line(:,1),s1_line(:,2));
P1 = InterX(s1_line2',full_branch1');
P1idx = find(min(P1(1,:)));
P1 = P1(:,P1idx);
if ismembertol(P1,top_branch1,1e-6)
    left_branch = top_branch1;
else
    left_branch = bot_branch1;
end    
idxp1s = dsearchn(s1_line2,P1');
idxp1d = dsearchn(left_branch,P1');

P2 = InterX(s1_line2',full_branch2');
P2idx = find(max(P2(1,:)));
P2 = P2(:,P1idx);
if ismembertol(P2,top_branch2,1e-6)
    right_branch = top_branch2;
else
    right_branch = bot_branch2;
end 
idxp2s = dsearchn(s1_line2,P2');
idxp2d = dsearchn(right_branch,P2');

dline = [left_branch(1:idxp1d-1,:);s1_line2(idxp1s+1:idxp2s-1,:);right_branch(idxp2d+1:end,:)];
z = [0,0,1];
[pt,dudt] = interparc(t2,dline(:,1),dline(:,2));
%% Identify which side of dline the xy points lie on
k = 1; m = 1;
tb = zeros(numel(xt));
top = zeros(ceil(0.75*numel(xt)),2);
bot = zeros(ceil(0.75*numel(xt)),2);
for i = 1:numel(xt)
    x0 = [xt(i),yt(i)];
    idx = dsearchn(pt,x0);
    nearest = pt(idx,:);
    tangent = [dudt(idx,:),0];
    pvec = [x0-nearest,0];
    p_cross_t = cross(pvec,tangent);
    if dot(p_cross_t,z) < 0
        top(k,:) = x0;
        k = k+1;
        tb(i) = -1;
    else
        bot(m,:) = x0;
        m = m+1;
        tb(i) = 1;
    end
end
%% Initialize
J = zeros(2,2,length(xt));
Sc = zeros(length(xt),2);
Y0 = cell(numel(xt),1);
if verify_plot
    [l_max,l_min,ev_max,ev_min] = CG_ev_field(vfield,tspan,x,y,T);    
    figure
    qquiver(X,Y,ev_min)
    hold on
    scatter2(degen_array,'filled');
    plot2(dline);
    scatter2(top); scatter2(bot);
    plot2(full_branch1); plot2(full_branch2);    
end

%% Perfrom coordinate transformation and obtain Jacobian
bs = polygon(2:3);
parfor i = 1:numel(xt)
    [s1_l,s2_l,J1,Y0{i}] = arclength_coord_tensorlines_tb(evfield2,intspan,X0(i,:),D,s1_arc,bs,tb(i));
    J(:,:,i) = J1;
    Sc(i,:) = [s1_l;s2_l];
end
%% Clean data and transform velocity for entire time epoch
fidx = find(abs(Sc(:,2))==0.001);
J(:,:,fidx) = NaN(2,2,numel(fidx));
Sc(fidx,:) = NaN(numel(fidx),2);
nvt = zeros(numel(xt),2,numel(tspan));
vt = zeros(size(nvt));
for i = 1:numel(xt)
    for k = 1:numel(tspan)
        vt = double_gyre(tspan(k),X0(i,:));
        nvt(i,:,k) = J(:,:,i)*vt';
        vt(i,:,k) = vt(i,:);        
    end
end