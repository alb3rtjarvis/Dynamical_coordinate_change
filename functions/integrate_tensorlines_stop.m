function [Y,s1,s2,error_flag] = integrate_tensorlines_stop(CG_eigvector_func,tspan,t0,y0,domain,T,Stop_line,bound_struct,init_direction,varargin)
% Nonadaptive integrator which computes tensorlines used to define new
% arclength coordinates from predifned grid in xy
% Because eigenvector fields are being used (e.g. e1 and -e1 would both 
% satisfy eigenvector equation) this integrator checks direction at each 
% step and forces continuity to obtain well defined tensorlines
%
%----------------------------------Input----------------------------------%
%
% - CG_eigvector_func - function representing eigenvector field which will 
% be interated in
%
% - tspan - vector containing t values at which integration will be
% performed
%
% - t0 - initial time
%
% - y0 - initial spatial point
%
% - domain - spatial domain over which integration is performed
%
% - T - finite integration time to compute Cauchy-Green strain tensor
%
% - Stop_line - arc obtained from orthogonal eigenvector field representing
% other arclength axis, used to stop integrator
%
% - bound_struct - structure containing boundaries of regions which contain
% degenerate points or other problematic strucure, arcs from this struct
% also stop integrator
%
% - init_direction - initial direction to begin integration
%
%--------------------------------Output-----------------------------------%
%
% - Y - array containing xy values of tensorline
%
% - s1,s2 - arclength coordinates determined by length of tensorline and
% arclength value on Stop_line
% 
% - error_flag - counter to determine if the integrator required a
% direction change and if error_flag exceeds 1, an error was made in the
% problem definition
%

xmin = domain(1,1);
xmax = domain(2,1);
ymin = domain(1,2);
ymax = domain(2,2);
if ~isnumeric(tspan)
  error('TSPAN should be a vector of integration steps.');
end

if ~isnumeric(y0)
  error('Y0 should be a vector of initial conditions.');
end

h = diff(tspan);
if any(sign(h(1))*h <= 0)
  error('Entries of TSPAN are not in order.') 
end  

try
  f0 = feval(CG_eigvector_func,tspan(1),y0,varargin{:});
catch
  msg = ['Unable to evaluate the ODEFUN at t0,y0. ',lasterr];
  error(msg);  
end  

y0 = y0(:);   % Make a column vector.
if ~isequal(size(y0),size(f0))
  error('Inconsistent sizes of Y0 and f(t0,y0).');
end

if ~exist('init_direction','var')
    init_direction = 0;
end

neq = length(y0);
N = length(tspan);
Y = zeros(neq,N);

% Method coefficients -- Butcher's tableau
%  
%   C | A
%   --+---
%     | B

C = [1/5; 3/10; 4/5; 8/9; 1];

A = [ 1/5,          0,           0,            0,         0
      3/40,         9/40,        0,            0,         0
      44/45        -56/15,       32/9,         0,         0
      19372/6561,  -25360/2187,  64448/6561,  -212/729,   0
      9017/3168,   -355/33,      46732/5247,   49/176,   -5103/18656];

B = [35/384, 0, 500/1113, 125/192, -2187/6784, 11/84];

% More convenient storage
A = A.'; 
B = B(:);      

nstages = length(B);
F = zeros(neq,nstages);
% direction  = zeros(1,nstages-1);
Y(:,1) = y0;
bline = [vertcat(bound_struct.x),vertcat(bound_struct.y)];      
stop_line = Stop_line(:,1:2);
line_tol = 5e-3;
bound_tol = 5e-3;
boundary_error = 1;
error_flag = 0;
edge_tol = 5e-3; %make smaller in future
while boundary_error >= 1
    if boundary_error > 2
        s1 = NaN;
        s2 = NaN;
        error_flag = 2;
        sprintf('Error - caught in direction loop, manual intervention required \n');
        break
    end
    for i = 2:N
      ti = tspan(i-1);
      hi = h(i-1);
      yi = Y(:,i-1);  
      % General explicit Runge-Kutta framework
      if i ~= 2 
        start_direction = sign(dot(F(:,end),feval(CG_eigvector_func,ti,yi,varargin{:})));
      else
        if boundary_error > 1 || init_direction == 1
            start_direction = -1;
        else
            start_direction = 1;
        end
      end
      F(:,1) = start_direction*feval(CG_eigvector_func,ti,yi,varargin{:});  
      for stage = 2:nstages
        tstage = ti + C(stage-1)*hi;
        ystage = yi + F(:,1:stage-1)*(hi*A(1:stage-1,stage-1));
        direction = sign(dot(F(:,stage-1),feval(CG_eigvector_func,tstage,ystage,varargin{:})));
        F(:,stage) = direction*feval(CG_eigvector_func,tstage,ystage,varargin{:});
      end  
      Y(:,i) = yi + F*(hi*B);
% Uncomment to display integration in real-time to check for errors in
% problem definition
% %       hold on
%       %plot(Y(1,1:i),Y(2,1:i))
%       scatter(Y(1,i),Y(2,i),7,'filled')
%       drawnow
% %       hold off
      bv = Y(:,i)' - bline;
      bound_dist = min(vecnorm(bv,2,2));
%     disp(bound_dist);     %%Uncomment to display distance to nearest boundary
      clear bline bv b_dist
      if bound_dist <= bound_tol || Y(1,i) >= xmax-edge_tol || Y(1,i) <= xmin+edge_tol || Y(2,i) >= ymax-edge_tol || Y(2,i) <= ymin+edge_tol
          boundary_error = boundary_error+1;
          error_flag = 1;
          clear bound_dist
          clear Y
          Y(:,1) = y0;
          break
      end
      stopv = Y(:,i)' - stop_line;
      [line_dist,idx] = min(vecnorm(stopv,2,2));
%       disp(line_dist); %%Uncomment to display distance to Stop_line
      if line_dist <= line_tol || Y(1,i) >= xmax || Y(1,i) <= xmin || Y(2,i) >= ymax || Y(2,i) <= ymin
          s1 = Stop_line(idx,3);
          s2 = tspan(i);
          boundary_error = 0;
          break
      end
    end
end
Y = Y.';
Y = Y(any(~Y==0,2),:);

end
