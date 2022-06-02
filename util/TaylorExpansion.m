% /*! @TaylorExpansion.m
% *************************************************************************
% <PRE>
% file.name       : TaylorExpansion.m
% related files   : TransFormSystem is LTI system transform
% function&ablity : TaylorExpansion and transform sys
% author          : gaodengwei
% version         : 1.00
% --------------------------------------------------------------------------------
% reference       : 
% --------------------------------------------------------------------------------
% record of modify :
% date          version     name         content
% 2016/02/28    1.00                     build
% </PRE>
% ********************************************************************************
%
% * right(c)
%
% *************************************************************************
% input:

% output: ploynomianl system:
%         zdothat£ºtransform to frame x-x0 with feedback
%         zdothatOrig:transform to frame x-x0 without feedback
%         xdothat£º feedback system 
%         xdothatOrig: open-loop system
%         system:
%         feedbacksys: closed-loop system
% *************************************************************************
function [zdothat, zdothatOrig, xdothat, xdothatOrig, feedbacksys]= TaylorExpansion(sys,order,varargin)

if nargin>2
    % get the new Cell ODE
    Gain = varargin{1}{1};
    u0 = varargin{1}{2};
else
    % if sys.K is set we use old matrix ODE
    Gain = polyniminalTrajectory(handle2ploy(sys.breaks,sys.K));
    % set the u0 to zeros and make it to ploynominal
    u0 = polyniminalTrajectory(handle2ploy(sys.breaks,@(t)zeros(sys.getNumInput,1)));
end
xtraj = sys.FunTraj;      % nominal trajectory
utraj = sys.FunContr;     % open loop control

% transfrom the coordinate frame: from state -> x-x0
ctf = PolynominalSystem(sys,[],[],[],[],eye(sys.getNumStates),-xtraj,0); % c sys
% x+x0 frame
dtf = PolynominalSystem(sys,[],[],[],[],eye(sys.getNumStates),xtraj,0);
 
ltvsys = PolynominalSystem(sys,[],[],[],[],Gain,u0,1); % control system with  (x-x0)

cInput = PolynominalSystem(sys,[],[],[],[],eye(sys.getNumInput),utraj,1);   % outputframe u+u0
%% closed-loop system: xdot = f(t,x,u0-K(x-x0))
% transfrom the coordinate frame: from x-x0 -> state
feedbacksys = cascade(ctf,ltvsys);              % state feedback system:y = 0*x+K*u-K*x0,x=0*x+0u+0
% transfrom the coordinate frame: from u+u0 -> state
feedbacksys = cascade(feedbacksys,cInput);      % output feedback system
%% x0 is trajectory
% dynamic computation
p_xu = [sys.p_x; sys.p_u]; % messploy

% taylor Approx
xdothat = @(t)build_poly(@feedbacksys.dynamics,t,xtraj.eval(t),u0.eval(t),order,p_xu);
xdothatOrig = @(t)build_poly(@sys.dynamics,t,xtraj.eval(t),utraj.eval(t),order,p_xu);
% trans z = x-x0 frame
zdothat = @(t) newdyn(t,sys,ctf,dtf,xdothat);
zdothatOrig = @(t) newdyn(t,sys,ctf,dtf,xdothatOrig);

end

function p=build_poly(fun,t,x0,u0,order,p_xu)
nX=length(x0);
nU=length(u0);
xu0=[x0;u0];
xu=TaylorVar.init(xu0,order);
x0=xu(1:nX);
u0=xu(nX+(1:nU));
p=getmsspoly(fun(t,x0,u0),p_xu-xu0);
end

% transform dynamic to new frame
function g = newdyn(t,sys,ctf,dtf,xdothat)
x = sys.p_x;
z = sys.p_x;
% p_t = msspoly('t',1);
[~,dc]=ctf.output(t,[],x);
d = dtf.output(t,[],z);
dcdt = dc(:,1); dcdx = dc(:,2:end);
g = dcdt + dcdx*xdothat(t);
g = subss(g,x,d);
end
