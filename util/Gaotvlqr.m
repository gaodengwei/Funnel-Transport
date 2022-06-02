function [K,V,B] = Gaotvlqr(obj,xtraj,utraj,Q,R,Vf,options)
% implements the time-varying linear quadratic regulator
% in the following xbar = x-x0, ubar = u-u0
typecheck(xtraj,'polyniminalTrajectory');
typecheck(utraj,'polyniminalTrajectory');

if (nargin<7), options=struct(); end 
 
nX = obj.getNumStates;
nU = obj.getNumInput;
if ~isfield(options,'tspan')
    tspan = obj.breaks;
else
    tspan = options.tspan;
end
if ~isfield(options,'N'), options.N=zeros(nX,nU); end

%  min_u(t) \int_0^T  x'Q{1}x + x'Q{2} + Q{3} + u'R{1}u + u'R{2} + R{3} + 2x'Nu
% subject to xdot = Ax + Bu + c

if isa(Q,'double')
    sizecheck(Q,[nX,nX]);
    Q = {ConstantTrajectory(Q),ConstantTrajectory(zeros(nX,1)),ConstantTrajectory(0)};
elseif isa(Q,'Trajectory')
    error('do not complete')
else
    error('do not complete')
end
if isa(R,'double')
    sizecheck(R,[nU,nU]);
    R = {ConstantTrajectory(R),ConstantTrajectory(zeros(nU,1)),ConstantTrajectory(0)};
elseif isa(R,'Trajectory')
    error('do not complete')
else
    error('do not complete')
end
if isa(Vf,'double')
    sizecheck(Vf,[nX,nX]);
    Qf = {Vf,zeros(nX,1),0};
elseif isa(Vf,'V_function')||isa(Vf,'quadraticV_function')
    Vf = inFrame(Vf,xtraj.eval(tspan(end))); 
    Qf{1}=Vf.S;
    Qf{2}=Vf.s1;
    Qf{3}=Vf.s2;
else
    error('can not find a lyapunov function form')
end
sizecheck(Qf,3);
typecheck(Qf{1},'double');
sizecheck(Qf{1},[nX,nX]);
typecheck(Qf{2},'double');
sizecheck(Qf{2},[nX,1]);
typecheck(Qf{3},'double');
sizecheck(Qf{3},1);

N=options.N;
if isa(N,'double')
    sizecheck(N,[nX,nU]);
    N = ConstantTrajectory(N);
elseif isa(N,'Trajectory')
    sizecheck(N,[nX,nU]);
else
    error('N must be a double or a trajectory');
end

xdottraj = fnder(xtraj);
% here use the sqrt method to solve time-varying lqr
Qf{1} = Qf{1}^(1/2);
S = cellODE(@(t,S)affineSdynamics(t,S,obj,Q,R,N,xtraj,utraj,xdottraj,options),tspan(end:-1:1),Qf);
S = flipToPT(S);
% back sqrt to nomal
S = recompS(S);
B = getBTrajectory(obj,S{1}.getBreaks(),xtraj,utraj,options);
% solve the feedback control
K = affineKsoln(S,R,B,N);

V = quadraticV_function(S{1},S{2},S{3});
V.x0 = xtraj;

end

% find the feedback control gain
function K = affineKsoln(S,R,B,N)
  Ri = inv(R{1});
  K{1} = -Ri*(B'*S{1}+N');
  K{2} = -.5*Ri*(B'*S{2} + R{2});
end
% make the B as a spline trajectory
function B = getBTrajectory(plant,ts,xtraj,utraj,options)
  nX = length(xtraj); nU = length(utraj);
  B = zeros([nX,nU,length(ts)]);
  for i=1:length(ts)
    x0=xtraj.eval(ts(i)); u0 = utraj.eval(ts(i));
  [~,df] = geval(@plant.dynamics,ts(i),x0,u0,options);    
    B(:,:,i) = df(:,nX+1+(1:nU));
  end
  B = polyniminalTrajectory(spline(ts,B));
end

function Sdot = affineSdynamics(t,S,plant,Qtraj,Rtraj,Ntraj,xtraj,utraj,xdottraj,options)
% solve the time varying lqr
x0 = xtraj.eval(t); u0 = utraj.eval(t); xdot0 = xdottraj.eval(t);

Q{1}=Qtraj{1}.eval(t); Q{2}=Qtraj{2}.eval(t); Q{3}=Qtraj{3}.eval(t);
R{1}=Rtraj{1}.eval(t); R{2}=Rtraj{2}.eval(t); R{3}=Rtraj{3}.eval(t);
Ri = inv(R{1});
N = Ntraj.eval(t);
nX = length(x0); nU = length(u0);

[xdot,df] = geval(@plant.dynamics,t,x0,u0,options);

A = df(:,1+(1:nX));
B = df(:,nX+1+(1:nU));
c = xdot - xdot0;

ss = inv(S{1}');
Sdot{1} = -.5*Q{1}*ss - A'*S{1} + .5*(N+S{1}*S{1}'*B)*Ri*(B'*S{1}+N'*ss);  % sqrt method  %#ok<MINV>
Sorig = S{1}*S{1}';

if (min(eig(Sorig))<0)
    warning('TVLQR: S is not positive definite');
end

rs = (R{2}+B'*S{2})/2;
Sdot{2} = -(Q{2} - 2*(N+Sorig*B)*Ri*rs + A'*S{2} + 2*Sorig*c);
Sdot{3} = -(Q{3}+R{3} - rs'*Ri*rs + S{2}'*c);
end

% cell ode
function output = cellODE(ODEFUN,tspan,A0)
if (~iscell(A0))
  error('A0 is not a cell array.  use normal ODE.');
end
if (numel(A0)~=length(A0)) 
  error('only handle cell vectors (nx1 or 1xn) so far');
end

y0=[];
for i=1:length(A0)
  y0 = [y0; reshape(A0{i},[],1)];
end

soln = ode45(@odefun,tspan,y0);
output = ODE45SolTrajectory(soln,cellfun(@size,A0,'UniformOutput',false));

function ydot = odefun(t,y)
    for j=1:length(A0)
      n = numel(A0{j});
      A{j} = reshape(y(1:n),size(A0{j}));
      y = y(n+1:end);
    end
    
    Adot = ODEFUN(t,A);

    ydot = [];
    for j=1:length(A0)
      ydot = [ydot; reshape(Adot{j},[],1)];
    end
end 
end

function S = recompS(Ssqrt)
  S{1} = Ssqrt{1}*Ssqrt{1}';
  S{2} = Ssqrt{2};
  S{3} = Ssqrt{3};
end