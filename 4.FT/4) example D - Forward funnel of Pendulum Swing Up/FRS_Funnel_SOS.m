function V = FRS_Funnel_SOS(sys,Vtraj0,ts,options) 

typecheck(Vtraj0,'quadraticV_function')
checkDependency('mosek');   
% checkDependency('sedumi');
N = length(ts);
x = sys.p_x;
u = sys.p_u;

if (~isfield(options,'rho0_tau')) options.rho0_tau=2; end
if (~isfield(options,'max_iterations')) options.max_iterations = 50; end
if (~isfield(options,'converged_tol')) options.converged_tol=.01; end
if (~isfield(options,'stability')) options.stability=false; end  % true implies that we want exponential stability
if (~isfield(options,'plot_rho')) options.plot_rho = true; end
if (~isfield(options,'lyap_parameterization')) options.lyap_parameterization = 'rho'; end
if (~isfield(options,'clean_tol')) options.clean_tol = 1e-6; end
if (~isfield(options,'rho0')) options.rho0 = .1*ones(length(ts),1); options.rho0(end) = 1; end
if (~isfield(options,'backoff_percent')) options.backoff_percent = 0; end 
if (~isfield(options,'degV')) options.degV = 2; end 
% Get the necessary variables

xdot = cell(1,N);
xdotOrig = cell(1,N);
V = cell(1,N);
Vdot = cell(1,N);
[Tayxdot, TayxdotOrig]= TaylorExpansion(sys,3,options.K); % taylor expantion in order

num_u = sys.getNumInput;
Vmin = zeros(N-1,1);
%% evaluate dynamics and Vtraj at every ts once (for efficiency/clarity)
for i=1:N
    V{i} = Vtraj0.getPoly(ts(i));
    xdot{i} = Tayxdot(ts(i));
    xdotOrig{i} = TayxdotOrig(ts(i));
    if (sys.getNumInput>0)   % zero all inputs
        xdot{i} = subs(xdot{i},sys.p_u,zeros(sys.getNumInput,1));
    end
    
    dVdt = Vtraj0.getPolyTimeDeriv(ts(i));
    Vdot{i} = diff(V{i},x)*xdot{i} + dVdt;
    Vdot{i} = clean(Vdot{i},options.clean_tol);
    V{i} = clean(V{i},options.clean_tol);
    
    % form new rho finding
    ui{i} = sys.FunContr.eval(ts(i)) + options.K{1}.eval(ts(i))*x;
    Vy{i} = Vtraj0.getPoly(ts(i));
    Vmin(i) =  minimumV(x,Vy{i});
    if i > 1
        Phi{i} = 0.5*eye(length(x));
    else
        Phi{i} = zeros(length(x),length(x));
    end
end
%%
% update the Vdot and V
%% Initialize rho with "tube" constraints
if (~isfield(options,'degL1'))
    options.degL1 = 4;%deg(Vdot{1},x);  % just a guess
end
utraj = sys.FunContr;
forig_u = xdotOrig;
%% Initialize rho
dts = diff(ts);
try
%     error
    load('rho.mat')
catch
    rho = rhoLineSearch(Vtraj0,Vy,utraj,ts,forig_u,Phi,dts,options,u,ui,x,Tayxdot);
end
rhodot = diff(rho)./dts;
% perform bilinear search the actual verification here

%% 
% rho = flip(rho);rhodot = diff(rho)./dts;
G = -.5*doubleSafe(subs(diff(diff(options.X0,x)',x),x,0*x)) / double(subs(options.X0,x,zeros(size(x))));
% Phase I: Find initial condition set
G_end = (Vtraj0.S.eval(ts(end)) + Phi{end})/rho(end);
G_end = G_end/5;

[rho_fac,L_tau] = computeRhoFac(Vy{1},G,rho(1),x);
rho_fac
L_tau
while rho_fac < 1.00
    % First step: Fix rho and V, search for L and u
    L1 = findL1(Vtraj0,Vy,rho,rhodot,utraj,ts,forig_u,Phi,options,u,ui,x);
%     plot(ts,rho);
%     drawnow;
%     hold on
    % Second step: Fix L, search for V and rho
    [Vy,rho,Phi] = findVRho_containG_fixed_end(G,Vtraj0,ts,forig_u,Phi,options,x,u,ui,L1,G_end,L_tau);
    rhodot = diff(rho)./dts;
 
    % Figure out what rho_fac is
    [rho_fac,L_tau] = computeRhoFac(Vy{1},G,rho(1),x);
    rho_fac
    L_tau
    % Plot stuff
    for k = 1:length(ts)
        S0(:,:,k) = eval(Vtraj0.S,ts(k))/rho(k);
        s10(:,k) = eval(Vtraj0.s1,ts(k))/rho(k);
        s20(k) = eval(Vtraj0.s2,ts(k))/rho(k);
        Phik(:,:,k) = double(Phi{k})/rho(k);
        S(:,:,k) = S0(:,:,k) + Phik(:,:,k);
    end
    STraj = polyniminalTrajectory(spline(ts(1:end),S));
    s1Traj = polyniminalTrajectory(spline(ts(1:end),s10));
    s2Traj = polyniminalTrajectory(spline(ts(1:end),s20));
    V = quadraticV_function(STraj,s1Traj,s2Traj); 
end
V = V*double(rho_fac);
options.rho0 = rho(1); 
%% Now, Phase II 
for k = 1:length(ts)
    S(:,:,k) = double(.5*subs(diff(diff(Vy{k},x)',x),x,0*x));
    s10(:,k) = double(subs(diff(Vy{k},x),x,0*x))';
    s20(k) = double(subs(Vy{k},x,0*x));  
end
STraj = polyniminalTrajectory(spline(ts(1:end),S));
s1Traj = polyniminalTrajectory(spline(ts(1:end),s10));
s2Traj = polyniminalTrajectory(spline(ts(1:end),s20));
V = quadraticV_function(STraj,s1Traj,s2Traj);  
% rho = ones(size(rho));
% rhodot = diff(rho)./dts;
% V = V/double(rho_fac(end));

B_vec = [];
options.rho0 = 1; 
options.degL1 = 4; 
options.degV = 2;
rho = initialRHO(V,ts,forig_u,dts,options,u,ui,x,rho);
rhodot = diff(rho)./dts;
for iter=1:options.max_iterations
    % First step: Fix rho and V, search for L and u
    [L1,Le] = findL(V,rho,rhodot,ts,forig_u,options,u,ui,x,G); 
    
    % Second step: Fix L, search for V and rho
    options.rho_all = rho;
    [V,rho] = findVRho_containG_vol(G,V,ts,forig_u,options,x,u,ui,L1,Le);
    rhodot = diff(rho)./dts;
    
    %% plot  
    tt = length(ts);
    for k = 1:tt
        S(:,:,k) = eval(V.S,ts(k))/rho(k);
        s10(:,k) = eval(V.s1,ts(k))/rho(k);
        s20(k) = eval(V.s2,ts(k))/rho(k);  
    end
    
    STraj = polyniminalTrajectory(spline(ts(1:end),S));
    s1Traj = polyniminalTrajectory(spline(ts(1:end),s10));
    s2Traj = polyniminalTrajectory(spline(ts(1:end),s20));
    V = quadraticV_function(STraj,s1Traj,s2Traj);
    
    figure
    options.plotdims = [1 2];
    options.inclusion = 'projection';
    options.x0 = sys.FunTraj;
      
    [~,B]= plot_myFunnel(sys,V,options);  
    disp([num2str(iter),'-region:',num2str(B)])
    title(['iteration ' num2str(iter)])
    drawnow;
    B_vec = [B_vec,B];
end
 

end

%%
function rho = initialRHO(Vtraj0,ts,forig_u,dts,options,u,ui,x,rho)

N = length(ts)-1;
% disp('Step 0: Initializing rho with bisection search...') 
 
try
for k = 1:N
    rhonow = rho(k);
    rhomin = 0.5*rhonow;
    rhomax = 2*rhonow; 
    rho(k+1) = fzero(@(rhonext) findRho(Vtraj0,rhonext,ts,forig_u,dts,options,u,ui,x,k,rhonow),[rhomin rhomax],optimset('TolX',1e-5));
    
%     rho(k+1) = 1.01*rho(k+1); % 1.01 To ensure feasibility
    disp([num2str(k),':  ',num2str(rho(k+1))])
end
catch
    keyboard
end

end 
function gamma = findRho(Vtraj0,rhonext,ts,forig_u,dts,options,u,ui,x,k,rho)
% Compute rhodot
rhodot = (rhonext - rho)/dts(k);

prog = spotsosprog;
prog = prog.withIndeterminate(x);
 
V0k = Vtraj0.getPoly(ts(k));
% V0dotk = Vtraj0.getPolyTimeDeriv(ts(k));
V0dotk = (Vtraj0.getPoly(ts(k+1))-V0k)/dts(k);

% Compute Vdot
Vdoty = diff(V0k,x)*subss(forig_u{k},u,ui{k}) + V0dotk;

% Clean stuff
V = clean(V0k,options.clean_tol);
Vdoty = clean(Vdoty,options.clean_tol);

% Declare multipliers
L1m = monomials(x,0:options.degL1);
[prog,L1] =  prog.newFreePoly(L1m); 
% Create gammas
[prog,gamma] = prog.newFree(1); 
% Declare SOS conditions
prog = prog.withSOS(-gamma*(x'*x)^(deg(Vdoty,x)/2) - Vdoty + rhodot - L1*(V-rho));
% prog = prog.withSOS(- Vdoty + rhodot + L1*(V-rho));
% Solve SOS program
pars = spot_sdp_default_options();
pars.verbose = 0;
sol = prog.minimize(-gamma,@spot_mosek,pars);

if strcmp(sol.info.solverInfo.itr.prosta,'PRIMAL_INFEASIBLE')
    gamma = -1.0;
else
    gamma = double(sol.eval(gamma));
    if strcmp(sol.info.solverInfo.itr.prosta,'UNKNOWN')
        gamma = -1.0;
    end
end
end

%%
function L1f = findL1(Vtraj0,Vy,rho,rhodot,utraj,ts,forig_u,Phi,options,u,ui,x)

N = length(ts)-1;
% disp('Step 1: Searching for multipliers and controller...')
 
% Optimized multipliers
L1f = cell(1,N);

for k = 1:N 
    prog = spotsosprog; % keyboard;
    prog = prog.withIndeterminate(x);
    
    Phidotk = (Phi{k+1} - Phi{k})/(ts(k+1) - ts(k));
    V0k = Vtraj0.getPoly(ts(k));
    V0dotk = Vtraj0.getPolyTimeDeriv(ts(k));
    
    % Compute Vdot
    Vdoty = diff(V0k,x)*subss(forig_u{k},u,ui{k}) + V0dotk + 2*x'*Phi{k}*subss(forig_u{k},u,ui{k}) + x'*Phidotk*x; 
    % Clean stuff
    V = clean(Vy{k},options.clean_tol);
    Vdoty = clean(Vdoty,options.clean_tol); 
    % Declare multipliers
    L1m = monomials(x,0:options.degL1);
    [prog,l1] = prog.newFree(length(L1m));
    L1 = l1'*L1m; 
    % Create gammas
    [prog,gamma] = prog.newPos(1); 
    % Declare SOS conditions
    prog = prog.withSOS(-gamma*(x'*x)^(deg(Vdoty,x)/2) - Vdoty + rhodot(k) + L1*(V-rho(k)));
    
    % Solve SOS program
    pars = spot_sdp_default_options();
    pars.verbose = 0;
    sol = prog.minimize(0,@spot_mosek,pars);
    
    % Optimized multipliers
    L1f{k} = sol.eval(L1);
end

end
function [rho_fac,L_tau] = computeRhoFac(V0,G,rho0,x)
prog = spotsosprog;

[prog,rho_fac] = prog.newPos(1);

[prog,L_tau] = prog.newPos(1);
prog = prog.withIndeterminate(x);
% Make sure r0 sublevel set of x'*G*x is contained in Phid{1} ellipse
prog = prog.withSOS(-L_tau*(V0 - rho0) - rho_fac + x'*G*x);

options = spot_sdp_default_options;
options.verbose = 0;
sol = prog.minimize(-rho_fac,@spot_mosek,options);

rho_fac = double(sol.eval(rho_fac));
L_tau = sol.eval(L_tau);
end

function [L1f,Lef]= findL(Vtraj0,rho,rhodot,ts,forig_u,options,u,ui,x,G) 
N = length(ts)-1; 
% Optimized multipliers
L1f = cell(1,N);
Lef = cell(1,N);
Lmonom = monomials(x,0:options.degL1);

% prog = spotsosprog;
% prog = prog.withIndeterminate(x); 
% [prog,L0] = prog.newFreePoly(Lmonom);
% prog = prog.withSOS(L0);
% prog = prog.withSOS(rho(1)-Vtraj0.getPoly(ts(1))-L0*(1-x'*G*x)); 

costfun = 0; 
for k = 1:N  
    prog = spotsosprog; % keyboard;
    prog = prog.withIndeterminate(x);
    [prog,L0] = prog.newFreePoly(Lmonom);
    V0k = Vtraj0.getPoly(ts(k));
    if k==1 
        prog = prog.withSOS(L0);
        prog = prog.withSOS(rho(1)-V0k-L0*(1-x'*G*x));
    end   
    V0dotk = (Vtraj0.getPoly(ts(k+1))-V0k)/(ts(k+1)-ts(k));
    % Compute Vdot
    Vdoty = diff(V0k,x)*subss(forig_u{k},u,ui{k}) + V0dotk; 
    % Clean stuff
    V = clean(V0k,options.clean_tol); 
    Vdoty = clean(Vdoty,options.clean_tol); 
    % Declare multipliers
    [prog,L1{k}] = prog.newFreePoly(Lmonom);  
    [prog,Le{k}] = prog.newFreePoly(Lmonom);  
    [prog,Sk] = prog.newPSD(length(x)); 
    prog = prog.withPSD(Sk); 
    % Create gammas 
    % Declare SOS conditions
    rhodot(k) = (rho(k+1)-rho(k))/(ts(k+1)-ts(k)); 
    % Create gammas
    [prog,gamma] = prog.newPos(1);   
%     gamma = 0;
    prog = prog.withSOS(-gamma*(x'*x)^(deg(Vdoty,x)/2) + rhodot(k) - Vdoty - L1{k}*(V-rho(k)));
    prog = prog.withSOS(1-x'*Sk*x-Le{k}*(rho(k)-V));
    prog = prog.withSOS(Le{k});
    % Solve SOS program
    pars = spot_sdp_default_options();
    pars.verbose = 0;
    [prog,tauk] = maxdet(prog,Sk);
%     costfun =  costfun-tauk;
    sol = prog.minimize(-tauk,@spot_mosek,pars); 
% %     % Optimized multipliers
    L1f{k} = sol.eval(L1{k});
    Lef{k} = sol.eval(Le{k});
end
% Solve SOS program
pars = spot_sdp_default_options();
pars.verbose = 0;
sol = prog.minimize(costfun/50,@spot_mosek,pars);
for k=1:N
    % Optimized multipliers
    L1f{k} = sol.eval(L1{k});
    Lef{k} = sol.eval(Le{k});
end

end

function [Vf,rho,Phi] = findVRho_containG_vol(G,Vtraj0,ts,forig_u,options,x,u,ui,L1,Le) 
N = length(ts)-1; 
% disp('Step 2: Searching for V and rho...') 
Vmonom = monomials(x,0:options.degV);
Lmonom = monomials(x,0:options.degL1);
prog = spotsosprog;
prog = prog.withIndeterminate(x);
% Declare rho
rho = msspoly(zeros(N+1,1)); 
   
% Make sure r0 sublevel set of x'*G*x is contained in Phid{1} ellipse  
[prog,V{1}] = prog.newFreePoly(Vmonom); 
[prog,L0] = prog.newFreePoly(Lmonom);
prog = prog.withSOS(L0);
prog = prog.withSOS(rho(1)-V{1}-L0*(1-x'*G*x)); 
[prog,V{N+1}] = prog.newFreePoly(Vmonom);

for k = N:-1:1 
     [prog,V{k}] = prog.newFreePoly(Vmonom); 
    V0dotk = (V{k+1}-V{k})/(ts(k+1) - ts(k));
     
    % Compute Vdot
    Vdoty = diff(V{k},x)*subss(forig_u{k},u,ui{k}) + V0dotk;
    
    % Clean stuff
    Vdoty = clean(Vdoty,options.clean_tol); 
    L1k = clean(L1{k},options.clean_tol);
    Lek = clean(Le{k},options.clean_tol);
    % Compute rhodot
    rhodot(k) = (rho(k+1)-rho(k))/(ts(k+1)-ts(k)); 
    % Declare SOS conditions
    prog = prog.withSOS(-Vdoty + rhodot(k) + L1k*(V{k}-rho(k))); 
    
    [prog,Sk] = prog.newPSD(length(x));
    prog = prog.withPSD(Sk);
    prog = prog.withSOS(1-x'*Sk*x-Lek*(rho(k)-V{k}));  
end

costfun = 0; 
for k = 1:(N+1) 
    if k == 1
        weight = 1; % 1
    else
        weight = 1; % 0.1;
    end 
    [prog,tauk] = maxdet(prog,Sk);
    costfun = costfun - weight*tauk;
end

costfun = costfun/(10); 
pars = spot_sdp_default_options();
pars.verbose = 0; 
sol = prog.minimize(costfun,@spot_mosek,pars);
 
if (options.backoff_percent > 0)
%     disp('Backing off now...')
    % Back off on objective
    costfun_opt = double(sol.eval(costfun));
    prog = prog.withPos(costfun_opt + (options.backoff_percent/100)*abs(costfun_opt) - costfun);
    sol = prog.minimize(0,@spot_sedumi,pars);
end

rho = double(sol.eval(rho));

for k = 1:(N+1)  
    Vf{k} = sol.eval(V{k});
end
end



function [V,rho,Phi] = findVRho_containG_fixed_end(G,Vtraj0,ts,forig_u,Phiold,options,x,u,ui,L1,G_end,L_tau)
N = length(ts)-1;
% disp('Step 2: Searching for V and rho...')
prog = spotsosprog;
prog = prog.withIndeterminate(x);
% Declare rho
[prog,rho] = prog.newPos(N+1);
rhodot = msspoly('r',N);
[prog,Phid{N+1}] = prog.newPSD(length(x)); 

% Normalization constraint
prog = prog.withEqs(trace(Phid{N+1}) - trace(Phiold{N+1}));

% Make V positive definite
[prog, Phislack] = prog.newPSD(length(x));
prog = prog.withEqs(Phislack - (Vtraj0.S.eval(ts(N+1)) + Phid{N+1}));

% Constrain V at end
V{N+1} = Vtraj0.getPoly(ts(N+1)) + x'*Phid{N+1}*x;
[prog,L0] = prog.newPos(1);
prog = prog.withSOS(V{N+1} - rho(end) - L0*(x'*G_end*x - 1));
[prog,Phid{1}] = prog.newPSD(length(x)); 

% Normalization constraint
prog = prog.withEqs(trace(Phid{1}) - trace(Phiold{1}));

% Make V positive definite
prog = prog.withPSD(Vtraj0.S.eval(ts(1)) + Phid{1});

% Make sure r0 sublevel set of x'*G*x is contained in Phid{1} ellipse
[prog,tau] = prog.newPos(1);

V{1} = Vtraj0.getPoly(ts(1)) + x'*Phid{1}*x;
prog = prog.withSOS(x'*G*x - tau - L_tau*(V{1} - rho(1)));

for k = N:-1:1
    if k > 1
        [prog,Phid{k}] = prog.newPSD(length(x));  
        % Normalization constraint
        prog = prog.withEqs(trace(Phid{k}) - trace(Phiold{k}));
    end 
    
    Phidotk = (Phid{k+1} - Phid{k})/(ts(k+1) - ts(k));
    V0k = Vtraj0.getPoly(ts(k));
    V0dotk = Vtraj0.getPolyTimeDeriv(ts(k));
    
    % Compute V
    V{k} = V0k + x'*Phid{k}*x;
    
    % Compute Vdot
    Vdoty = diff(V0k,x)*subss(forig_u{k},u,ui{k}) + V0dotk + 2*x'*Phid{k}*subss(forig_u{k},u,ui{k}) + x'*Phidotk*x;
    
    % Clean stuff
    Vdoty = clean(Vdoty,options.clean_tol);
    L1k = clean(L1{k},options.clean_tol);
    
    % Compute rhodot
    rhodot(k) = (rho(k+1)-rho(k))/(ts(k+1)-ts(k));
    
    % Declare SOS conditions
    prog = prog.withSOS(-Vdoty + rhodot(k) + L1k*(V{k}-rho(k)));
    
    % Make V positive definite
    [prog,Phislack] = prog.newPSD(length(x));
    prog = prog.withEqs(Phislack - (Vtraj0.S.eval(ts(k)) + Phid{k}));
    
end

costfun = -tau;

pars = spot_sdp_default_options();
pars.verbose = 0;

sol = prog.minimize(costfun,@spot_mosek,pars); 
rho = double(sol.eval(rho));

for k = 1:(N+1)
    Phi{k} = double(sol.eval(Phid{k}));
    V0k = Vtraj0.getPoly(ts(k));
    V{k} = V0k + x'*Phi{k}*x;
end
end



function y=doubleSafe(x)
y=double(x);
if (~isa(y,'double')) error('double failed'); end
end

function rho = rhoLineSearch(Vtraj0,Vy,utraj,ts,forig_u,Phi,dts,options,u,ui,x,psys)

N = length(ts)-1;
% disp('Step 0: Initializing rho with bisection search...')

rho = zeros(N+1,1);
rho(1) = options.rho0(1);

for k = 1:N
    rhonow = rho(k);
    rhomin = 0.01*rhonow;
    rhomax = 100*rhonow;
    rho(k+1) = fzero(@(rhonext) checkRho(Vtraj0,Vy,rhonext,utraj,ts,forig_u,Phi,dts,options,u,ui,x,k,rhonow,psys),[rhomin rhomax],optimset('TolX',1e-5));
    rho(k+1) = 1.01*rho(k+1); % 1.01 To ensure feasibility
    disp([num2str(k),':  ',num2str(rho(k+1))])
end

end

function gamma = checkRho(Vtraj0,Vy,rhonext,utraj,ts,forig_u,Phi,dts,options,u,ui,x,k,rho,psys)
% Compute rhodot
rhodot = (rhonext - rho)/dts(k);

prog = spotsosprog;
prog = prog.withIndeterminate(x);

Phidotk = (Phi{k+1} - Phi{k})/(ts(k+1) - ts(k));
V0k = Vtraj0.getPoly(ts(k));
V0dotk = Vtraj0.getPolyTimeDeriv(ts(k));

% Compute Vdot
Vdoty = diff(V0k,x)*subss(forig_u{k},u,ui{k}) + V0dotk + 2*x'*Phi{k}*subss(forig_u{k},u,ui{k}) + x'*Phidotk*x;

% Clean stuff
V = clean(Vy{k},options.clean_tol);
Vdoty = clean(Vdoty,options.clean_tol);

% Declare multipliers
L1m = monomials(x,0:options.degL1);
[prog,l1] = prog.newFree(length(L1m));
L1 = l1'*L1m;

% Create gammas
[prog,gamma] = prog.newFree(1);

% Declare SOS conditions
prog = prog.withSOS(-gamma*(x'*x)^(deg(Vdoty,x)/2) - Vdoty + rhodot + L1*(V-rho));

% Solve SOS program
pars = spot_sdp_default_options();
pars.verbose = 0;
sol = prog.minimize(-gamma,@spot_mosek,pars);

if strcmp(sol.info.solverInfo.itr.prosta,'PRIMAL_INFEASIBLE')
    gamma = -1.0;
else
    gamma = double(sol.eval(gamma));
    if strcmp(sol.info.solverInfo.itr.prosta,'UNKNOWN')
        gamma = -1.0;
    end
end
end


function [Vmin,b] = minimumV(x,V)
if (deg(V,x)>2)
    prog = mssprog;
    [prog,slack] = new(prog,1,'free');
    prog.sos = slack + V;
    [prog,info] = sedumi(prog,slack,0);
    Vmin = -doubleSafe(prog(slack));
else
    H = doubleSafe(0.5*diff(diff(V,x)',x));
    b = -0.5*(H\doubleSafe(subs(diff(V,x),x,0*x)'));
    Vmin = subs(V,x,b);
end
end