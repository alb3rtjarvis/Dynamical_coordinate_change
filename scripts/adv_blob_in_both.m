%% Compute solutions in xy
Opt = odeset('RelTol',1e-9,'AbsTol',1e-10);
addpath('../functions')
for k = 1:length(X0)
    [~,y_sol{k}] = ode45(vfield,[0,T/2],X0(k,:),Opt);
    disp(k);
end
%% Identify which trajectories stay in each piece of our domain
for k = 1:numel(y_sol) 
    for i = 1:numel(polygon)
        in_y{i,k} = inpolygon(y_sol{k}(:,1),y_sol{k}(:,2),polygon(i).x,polygon(i).y);
    end
end

%% Remove the points which fall in the degenerate domains
for k = 1:numel(y_sol)
    for i = 1:numel(polygon)-1
        in_y{1,k}(in_y{i+1,k}) = 0;
    end
end
%% Identify the index of trajectories which lie completely in our domain
i = 1;
for k = 1:numel(y_sol)
    if all(in_y{1,k})
        stay_idx(i) = k;
        i = i+1;
        disp(k);
    end
end

%% Identify indicies of boundary of boxes in xy to be advected
b(1,:) = [0.68,0.76];
b(2,:) = [0.52,0.78];
b(3,:) = [0.36,0.74];
b(4,:) = [1.7,0.74];
delta = 0.02;
for k = 1:length(b)
    box_xy(:,:,k) = aux_grid8(b(k,:),delta);
    for j = 2:size(box_xy,1)
        b_box_idx(j-1,k) = find(ismembertol(X0,box_xy(j,:,k),1e-10,'ByRows',true));
    end
    box_s(:,:,k) = Sc(b_box_idx(:,k),:);
end
box_xy(1,:,:) = [];
%% Identify points in xy and s1s2 that lie in initial boxes to be advected
xx = linspace(0,2,901);
yy = linspace(0,1,451);
[XX,YY] = meshgrid(xx,yy);
xxq = XX(:); yyq = YY(:);
s1min = min(Sc(:,1)); s1max = max(Sc(:,1));
s2min = min(Sc(:,2)); s2max = max(Sc(:,2));
ss1 = linspace(s1min,s1max,901);
ss2 = linspace(s2min,s2max,451);
[S1,S2] = meshgrid(ss1,ss2);
s1q = S1(:); s2q = S2(:);
X0b = cell(length(b),1);
for k = 1:length(b)
    in_b = inpolygon(xxq,yyq,box_xy(:,1,k),box_xy(:,2,k));
    X0b{k} = [xxq(in_b),yyq(in_b)];
    in_bs = inpolygon(s1q,s2q,box_s(:,1,k),box_s(:,2,k));    
    S0b{k} = [s1q(in_bs),s2q(in_bs)];
    clear in_b in_bs
end

%% Advect initial boxes in xy and s1s2
nt = 301;
for j = 1:size(X0b,2)
    for k = 1:size(X0b{j},1)
    [~,y_sol_bt(:,:,k)] = ode45(vfield,linspace(0,T/2,nt),X0b{j}(k,:),Opt);
    end
    y_sol_b{j} = y_sol_bt;
    clear y_sol_bt;
    disp(j);
end

for j = 1:size(S0b,2)
    for k = 1:size(S0b{j},1)
    [~,S_sol_bt(:,:,k)] = ode45(S_interp,linspace(0,T/2,nt),S0b{j}(k,:),Opt);
    end
    S_sol_b{j} = S_sol_bt;
    clear S_sol_bt;
    disp(j);
end
%% Reindex for simpler plotting
% for j = 1:size(X0b,3)
%     for k = 1:size(y_sol_c,1)
%         c_adv(:,:,k,j) = y_sol_c(k,:,:,j);
%     end
% end

% for j = 1:size(X0b_orig,3)
%     for k = 1:size(y_sol_c_orig,1)
%         c_adv_orig(:,:,k,j) = y_sol_c_orig(k,:,:,j);
%     end
% end
%% Create succesive figures of advected circles in xy for movie
sz = 36;
h = figure('units','pixels','position',[0 0 1920 1080]);
for k = 1:nt
    plot2(X_out1f)
    hold on
    plot2(X_out1b)
    plot(x_left,y,x_right,y)
    plot(x,y_bot,x,y_top)
    plot(d1x_line,d1y_line_top,d1x_line,d1y_line_bot);
    plot(d2x_line,d2y_line_top,d2x_line,d2y_line_bot);
    plot(d1vx,d1vy,d2vx,d2vy)
    
    scatter(y_sol_b{1}(k,1,:),y_sol_b{1}(k,2,:),'b','filled') 
    scatter(y_sol_b{2}(k,1,:),y_sol_b{2}(k,2,:),'r','filled') 
    scatter(y_sol_b{3}(k,1,:),y_sol_b{3}(k,2,:),'g','filled') 
    scatter(y_sol_b{4}(k,1,:),y_sol_b{4}(k,2,:),'c','filled') 
    drawnow
    mov(k) = getframe(gcf);
    disp(k);
	hold off
end

%% Create succesive figures of advected circles in xy for movie
sz = 36;
%load('transform_gridlines.mat');
% h = figure('units','pixels','position',[0 0 1920 1080]);
% hold on
% for j = 1:numel(syline)
%     plot2(syline{j},'r')
% end
% 
% for j = 1:numel(sxline)
%     plot2(sxline{j},'b')
% end
for k = 1:nt
h = figure('units','pixels','position',[0 0 1920 1080]);
hold on
    for j = 1:numel(syline)
        plot2(syline{j},'r')
    end

    for j = 1:numel(sxline)
        plot2(sxline{j},'b')
    end    
    scatter(S_sol_b{1}(k,1,:),S_sol_b{1}(k,2,:),'b','filled') 
    scatter(S_sol_b{2}(k,1,:),S_sol_b{2}(k,2,:),'r','filled') 
    scatter(S_sol_b{3}(k,1,:),S_sol_b{3}(k,2,:),'g','filled') 
    scatter(S_sol_b{4}(k,1,:),S_sol_b{4}(k,2,:),'c','filled') 
    drawnow
    mov(k) = getframe(gcf);
    disp(k);
    close all
end
%% Make movie
fps = 50;
obj = VideoWriter('dg_s_vel','MPEG-4');
obj.FrameRate = fps;
open(obj);
writeVideo(obj,mov);
close(obj);

%% Clean and reorganize data to create interpolant of s1s2 velocity 
nnan_idx = ~isnan(Sc(:,1));
slen = length(nnan_idx(nnan_idx));
tlen = length(tspan);
P = zeros(slen*tlen,3);
for k = 1:tlen
    P(1+(k-1)*slen:k*slen,:) = [Sc(nnan_idx,1),Sc(nnan_idx,2),tspan(k)*ones(slen,1)];
    NV1(1+(k-1)*slen:k*slen,1) = nvt(nnan_idx,1,k);
    NV2(1+(k-1)*slen:k*slen,1) = nvt(nnan_idx,2,k);
end
%% Create interpolant of s1s2 velocity
su_interp = scatteredInterpolant(P,NV1);
sv_interp = scatteredInterpolant(P,NV2);
% tvel = @tv_interp(Sc(:,1),Sc(:,2),tspan,nvt)
%% Advect corresponding initial blobs in s1s2
S_interp = @(t,x) tv_interp(t,x,su_interp,sv_interp);
for j = 1:size(S0b,3)
    for k = 1:size(S0b,1)
    [~,s_circ_sol(:,:,k,j)] = ode45(S_interp,linspace(0,T/2,301),S0b(k,:,j),Opt);
    end
    disp(j);
end
%% Reindex for simpler plotting
for j = 1:size(S0b,3)
    for k = 1:size(s_circ_sol,1)
        s_adv(:,:,k,j) = s_circ_sol(k,:,:,j);
    end
end
%% Create succesive figures of advected blobs in s1s2 for movie
hs = figure('units','pixels','position',[0 0 1920 1080]);
for k = 1:size(s_adv,3)
    scatter2(Sc);
    hold on
    scatter2(s_adv(:,:,k,1),sz,'b','filled') 
    scatter2(s_adv(:,:,k,2),sz,'r','filled') 
    scatter2(s_adv(:,:,k,3),sz,'g','filled') 
    scatter2(s_adv(:,:,k,4),sz,'c','filled') 
    drawnow
    %mov(k) = getframe(gcf);
    disp(k);
	hold off
end

%% Functions
function transformed_vel_interp = tv_interp(t_in,S_in,su_interp,sv_interp)
    s1 = S_in(1);
    s2 = S_in(2);
    su = su_interp(s1,s2,t_in);
    sv = sv_interp(s1,s2,t_in);
    transformed_vel_interp = [su;sv];
end