function [V,rho] = regionOfAttraction(sys,varargin)
% Estimates the region of attraction.

ok_sedumi = checkDependency('sedumi');
ok_mosek = checkDependency('mosek');
% ok_mosek = 0;
% ok_yalmip = checkDependency('yalmip');
if ~ok_sedumi && ~ok_mosek
    error('You need either MOSEK or SeDuMi installed to use this function.');
end

%% get Lyapunov candidate
num_x = sys.getNumStates;
num_u = sys.getNumInput;
t = msspoly('t');
x = msspoly('x',num_x);
u = msspoly('u',num_u);

% lyapunov equation
V = varargin{1};

% ploynominal or not
try
    f= sys.dynamics(t,x,u);
catch
    f = TaylorApproxSys(sys,sys.xG,sys.uG,3);
end
% default is origin(if not trans sys to V frame)
if ~all(V.x0==0)
    f = TransSystem(sys,V);
end
%% zero all inputs
if (sys.getNumInput>0)
    f = subs(f,u,zeros(sys.getNumInput,1));
end

%% handle options
if (nargin>2) options=varargin{2};
else options=struct(); end
if (~isfield(options,'method'))
    options.method={'levelset'};% bilinear 
elseif (~iscell(options.method))
    options.method={options.method};
end
if (~isfield(options,'degV')) options.degV = 4; end
if (~isfield(options,'max_iterations')) options.max_iterations=10; end
if (~isfield(options,'converged_tol')) options.converged_tol=.01; end
if (~isfield(options,'optimize')) options.optimize=true; end
if (~isfield(options,'clean_tol')) options.clean_tol = 1e-6; end % tolerance for cleaning small terms
if (~isfield(options,'numSamples')) options.numSamples = 10^(num_x+1); end

if (~isfield(options,'degL1'))
    options.degL1 = options.degV-1 + deg(f,x);  % just a guess
end
if (~isfield(options,'degL2'))
    options.degL2 = options.degL1;
end

% Choose sdp solver
if ok_mosek
    options.solver = @spot_mosek;
else
    options.solver = @spot_sedumi;
end

for i=1:length(options.method)
    %% compute level set
    switch (lower(options.method{i}))
        case 'bilinear'
            [V,rho] = bilinear(V,f,options);
        case 'levelset'
            [V,rho] = levelSetMethod(V,f,options);
        case 'stochasticsearch'
            f = sys.dynamics(t,x,u);
            V = stochasticsearch(V,f,u,options); 
        case 'parsearch'
            f = sys.dynamics(t,x,u);
%             [V,rho] = parsearch(V,f,u,options); % sample based ROA
            [V,rho] = parsearchCT(V,f,u,options); % CT based ROA
        case 'bilinearequality'
            V = bilinearequality(V,f,options);
        otherwise
            error(['don''t know method: ', options.method]);
    end
end

end

function [T,Vbal,fbal,S,A] = balance(x,V,f,S,A)
if (nargin<4 || isempty(S))
    S=.5*doubleSafe(subs(diff(diff(V,x)',x),x,0*x));  % extract Hessian
end
if (nargin<5 || isempty(A))
    A = doubleSafe(subs(diff(f,x),x,0*x));
end

[T,D] = balanceQuadForm(S,(S*A+A'*S));

if (nargout>1)
    Vbal=subs(V,x,T*x);
    if (nargout>2)
        fbal=inv(T)*subs(f,x,T*x);
    end
end
end
 
%% for the bilinear search
function [V,rho] = bilinear(V0,f,options)

x = V0.x; 
V=V0; 
[T,V0bal,fbal,S0,A] = balance(x,V0.getPoly,f);
rho = 1;

L1monom = monomials(x,0:options.degL1);
L2monom = monomials(x,0:options.degL2);
Vmonom = monomials(x,0:options.degV);

vol=0;
Vpoly = V.getPoly;
for iter=1:options.max_iterations 
    last_vol = vol; 
    % balance on every iteration (since V and Vdot are changing):
    [T,Vbal,fbal]=balance(x,Vpoly,f,S0/rho,A);
    V0bal=subs(V0.getPoly,x,T*x);
    
    [L1,sigma1] = findL1(x,fbal,Vbal,L1monom,options);
    L2 = findL2(x,Vbal,V0bal,rho,L2monom,options);
    [Vbal,rho] = optimizeV(x,fbal,L1,L2,V0bal,sigma1,Vmonom,options);
    vol = rho;
    Vpoly = clean(subs(Vbal,x,inv(T)*x));
    disp([num2str(iter),':',num2str(rho)]);
    
%     y2 = getLevelSet(x,Vpoly,options);
%     plot3(y2(1,:),y2(2,:),repmat(3,1,size(y2,2)),'r','LineStyle','-','LineWidth',2);

    V = Vpoly;
    % check for convergence
    if ((vol - last_vol) < options.converged_tol*last_vol)
        break;
    end
end 
end

function [L1,sigma1] = findL1(x,f,V,Lxmonom,options)
prog = spotsosprog;
prog = prog.withIndeterminate(x);

% construct multipliers for Vdot
[prog,L1] = prog.newFreePoly(Lxmonom);

% construct Vdot
Vdot = clean(diff(V,x)*f);

% construct slack var
[prog,sigma1] = prog.newPos(1);

% setup SOS constraints
prog = prog.withSOS(-Vdot + L1*(V - 1) - sigma1*V);
prog = prog.withSOS(L1);

% run SeDuMi/MOSEK and check output
solver = options.solver;
options = spot_sdp_default_options();
sol = prog.minimize(-sigma1,solver,options);

if sol.status == spotsolstatus.STATUS_SOLVER_ERROR
    error('The solver threw an internal error.');
end
if ~sol.isPrimalFeasible
    error('Problem looks primal infeasible');
end

if ~sol.isDualFeasible
    error('Problem looks dual infeasible. It is probably unbounded. ');
end

L1 = clean(sol.eval(L1));
sigma1 = sol.eval(sigma1);
end

function L2 = findL2(x,V,V0,rho,Lxmonom,options)
prog = spotsosprog;
prog = prog.withIndeterminate(x);

% construct multipliers
[prog,L2] = prog.newFreePoly(Lxmonom);

[prog,slack] = prog.newPos(1);

prog = prog.withSOS(-(V-1) + L2*(V0-rho));
prog = prog.withSOS(L2);

solver = options.solver;
options = spot_sdp_default_options();
options.verbose = 0;
sol = prog.minimize(slack,solver,options);% keyboard;

if ~sol.isPrimalFeasible
    error('Problem looks primal infeasible');
end

if ~sol.isDualFeasible
    error('Problem looks dual infeasible. It is probably unbounded. ');
end

L2 = clean(sol.eval(L2));
end

function [V,rho]=optimizeV(x,f,L1,L2,V0,sigma1,Vxmonom,options)
prog = spotsosprog;
prog = prog.withIndeterminate(x);

% construct V
[prog,V] = prog.newFreePoly(Vxmonom);
Vdot = diff(V,x)*f;

% construct rho
[prog,rho] = prog.newPos(1);

% setup SOS constraints
prog = prog.withSOS(-Vdot + L1*(V - 1) - sigma1*V/2);
prog = prog.withSOS(-(V-1) + L2*(V0 - rho));
prog = prog.withSOS(V);

% run SeDuMi/MOSEK and check output
solver = options.solver;
options = spot_sdp_default_options();
options.verbose=0;
sol = prog.minimize(-rho,solver,options);

if ~sol.isPrimalFeasible
    error('Problem looks primal infeasible');
end

if ~sol.isDualFeasible
    error('Problem looks dual infeasible. It is probably unbounded.');
end

V = sol.eval(V);
rho = double(sol.eval(rho));
end


%% Pablo's method (jointly convex in rho and lagrange multipliers)
function [V,rho] = levelSetMethod(V0,f,options)

x = V0.x;
[T,V,f] = balance(x,V0.getPoly,f);

%% compute Vdot
Vdot = diff(V,x)*f;

V = clean(V);
Vdot = clean(Vdot);

% check Hessian Vdot at origin, to make sure it's negative def.
H=.5*doubleSafe(subs(diff(diff(Vdot,x)',x),x,0*x));  % extract Hessian
if (~isPositiveDefinite(-H)) error('Vdot must be negative definite at the origin'); end

prog = spotsosprog;
prog = prog.withIndeterminate(x);
Lmonom = monomials(x,0:options.degL1);
%  Lmonom = hermite_basis(monomials(x,0:options.degL1));

[prog,rho] = prog.newFree(1);

[prog,L] = prog.newFreePoly(Lmonom);

prog = prog.withSOS((x'*x)^floor((options.degL1 + deg(Vdot)-deg(V))/2)*(V - rho) +  L*Vdot);

solver = options.solver;
options = spot_sdp_default_options();
options.verbose = 0;
sol = prog.minimize(-rho,solver,options);

if ~sol.isPrimalFeasible
    error('Problem looks primal infeasible');
end

if ~sol.isDualFeasible
    error('Problem looks dual infeasible. It is probably unbounded. ');
end

rho = doubleSafe(sol.eval(rho));
if (rho<=0) error('optimization failed'); end

V = V/rho;

%% undo balancing
V = subs(V,x,inv(T)*x); 


end

%% uncertain search
function V = stochasticsearch(V0,f,u,options)
x = V0.x; 
V = V0;  

gW = [u-options.box(:,1);options.box(:,2)-u];
degLw = options.degV -1 + deg(gW,u);  % just a guess

L1monom = monomials([x;u],0:options.degL1); 
Lwmonom = [monomials([x;u],0:degLw)]; 
L2monom = monomials(x,0:options.degL2);
Vmonom = monomials(x,0:options.degV);
rho = 1;
vol=0;
V = V.getPoly; 
p = 2*V;
for iter=1:options.max_iterations
    last_vol = vol; 
    p = 1.01*V; % backoff
    % balance on every iteration (since V and Vdot are changing):   
    [L1,Lw] = findLw(x,u,f,V,L1monom,Lwmonom,gW,options);
    L2 = findLw2(x,V,p,1,L2monom,options);
    [V,rho] = optimizeVw(x,u,f,L1,L2,Lw,Vmonom,p,gW,options);   
     iter 
     vol = rho;
    % check for convergence
    if (abs(vol - last_vol) < 0.01*options.converged_tol*last_vol)
        break;
    end
end
end

function [L1,L2] = findLw(x,w,f,V,Lxmonom,Lwmonom,gW,options)  
prog = spotsosprog;
prog = prog.withIndeterminate([x;w]);

% construct multipliers
[prog,L1] = prog.newFreePoly(Lxmonom);
Lw = [];
for i=1:length(gW)
    [prog,Lwmom] = prog.newFreePoly(Lwmonom); 
    Lw = [Lw;Lwmom];
end 
V = clean(V,options.clean_tol);
Vdot = clean(diff(V,x)*f,options.clean_tol); 

[prog,gamma0] = prog.newPos(1); 
prog = prog.withSOS( -gamma0*(x'*x)^2- Vdot - L1*(1- V) - Lw'*gW);%
prog = prog.withSOS(Lw);
prog = prog.withSOS(L1);
 
solver = options.solver;
pars = spot_sdp_default_options();
 
sol = prog.minimize(-gamma0,solver,pars);

if sol.status == spotsolstatus.STATUS_SOLVER_ERROR
    error('The solver threw an internal error.');
end
if ~sol.isPrimalFeasible
    error('Problem looks primal infeasible');
end

if ~sol.isDualFeasible
    error('Problem looks dual infeasible. It is probably unbounded. ');
end

L1 = sol.eval(L1);
L2 = sol.eval(Lw); 

end

function L2= findLw2(x,V,p,rho,L2monom,options) 

prog = spotsosprog;
prog = prog.withIndeterminate(x); 
 
[prog,L2] = prog.newFreePoly(L2monom);
  
prog = prog.withSOS(1-V-L2*(rho-p));
prog = prog.withSOS(L2); 

solver = options.solver;
options = spot_sdp_default_options();

sol = prog.minimize(0,solver,options);

if ~sol.isPrimalFeasible
    error('Problem looks primal infeasible');
end

if ~sol.isDualFeasible
    error('Problem looks dual infeasible. It is probably unbounded. ');
end

L2 = sol.eval(L2);
end

function [V,rho]=optimizeVw(x,w,f,L1,L2,Lwmonom,Vxmonom,p,gW,options)
 
prog = spotsosprog;
prog = prog.withIndeterminate([x;w]);

% Create Phi and new V 
[prog,V] = prog.newFreePoly(Vxmonom); 
Vdot = diff(V,x)*f;  
 
% construct uncertainty
 Lw = [];
for i=1:length(gW)
    [prog,Lwmom] = prog.newFreePoly(Lwmonom); 
    Lw = [Lw;Lwmom];
end    
  
[prog,rho] = prog.newPos(1);  

% setup SOS constraints
prog = prog.withSOS(- Vdot - L1*(1- V) - Lw'*gW);
prog = prog.withSOS(1-V - L2*(rho-p));
prog = prog.withSOS(V);
prog = prog.withSOS(Lw); 
% run SeDuMi/MOSEK and check output
 
solver = options.solver;
options = spot_sdp_default_options();
sol = prog.minimize(-rho,solver,options);

if ~sol.isPrimalFeasible
    error('Problem looks primal infeasible');
end

if ~sol.isDualFeasible
    error('Problem looks dual infeasible. It is probably unbounded. ');
end

V = clean(sol.eval(V));
rho = double(sol.eval(rho));
end

 
function y=doubleSafe(x)
y=double(x);
if (~isa(y,'double')) error('double failed'); end
end
