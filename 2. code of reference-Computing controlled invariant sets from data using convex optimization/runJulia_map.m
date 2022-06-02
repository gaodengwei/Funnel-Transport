% Code for the paper "Data-driven convex computation of invariant sets for nonlinear dynamical systems"
% Example: Julia map
% Author: Milan Korda, 2019

% Readme: Needs a linear programming solver installed. The numerical
% examples in the paper were run with Gurobi, but other solvers will work
% as well (mosek, cplex, linprog etc)

clear all
close all
checkDependency('yalmip') 
checkDependency('mosek')
addpath(genpath('./Resources'));

rng(3141414)

% Linear programming solver
LP_SOLVER = 'linprog'; % Change to another LP solver if you have it, it's probably better

% Dynamics
cc = [-0.7 ; 0.2];
f = @(x)([ x(1,:).^2 - x(2,:).^2 + cc(1) ;
           2*x(1,:).*x(2,:) + cc(2) ]); % Julia recursion
n = 2;

% Discount factor
alpha = 0.6;

% Basis
basis = 'monomial';
switch basis
    case 'monomial'
        d = 10;
        basisEval = @(x)(evalMonBasis(d,x));
        n_basis = numMonTot(n,d);
    case 'rbf'
        n_basis = 400; 
        cent = rand_sphere(n,n_basis);           
        basisEval = @(x)( rbf(x,cent,'thinplate') );   
end

% Generate data
n_samp = 3e4;

CONSERVATIVE = 0;
if(CONSERVATIVE == 0)
    % Everything for optimization
    x_samp = rand_sphere(n,n_samp);
    x_samp_plus = f(x_samp);
else
    % Half and half for optimization and validation
    x_samp = rand_sphere(n,round(n_samp / 2));
    x_samp_valid = rand_sphere(n,round(n_samp / 2));
    x_samp_plus = f(x_samp);                
    x_samp_valid_plus = f(x_samp_valid);    
end


% Projection function
proj_X = @(x)projBall(x);

% x^+
x_samp_plus_proj = proj_X(x_samp_plus);

% \bar l
l_samp = min(sqrt(sum( (x_samp_plus_proj - x_samp_plus).^2 )),1);

% Evaluation of the basis on samples
basis_samp = basisEval(x_samp);
basis_samp_plus = basisEval(x_samp_plus_proj);



%% Set up and solve optimization problem

% Decision variable
cv = sdpvar(n_basis,1); 

v_samp = cv'*basis_samp;
v_samp_plus = cv'*basis_samp_plus;

% Constraints
con = [v_samp <= l_samp + alpha*v_samp_plus]; % v <= l + alpha v\circ f on data

% Objective
n_samp_obj = 1e5;
x_samp_obj = rand_sphere(n,n_samp_obj); % Samples for monte carlo approximation of the objective = - int_X v(x) dx
v_samp_obj = cv'*basisEval(x_samp_obj);
obj = -sum(v_samp_obj) / n_samp_obj;

% Solver parameters
switch lower(LP_SOLVER)
    case 'gurobi'
        prec = 1e-6;
        options = sdpsettings('solver','gurobi','gurobi.Method',2, 'gurobi.Crossover', 0,...
            'gurobi.BarConvTol', prec, 'gurobi.FeasibilityTol',prec, 'gurobi.OptimalityTol',prec);
    otherwise
        options = sdpsettings('solver', LP_SOLVER);
end

% Solve
diag = optimize(con,obj,options);
if(diag.problem ~= 0)
    warning(diag.info)
end

%% Extract solution

cv = double(cv);

fprintf('Objective value = %f \n', -double(obj))
v = @(x)(cv'*basisEval(x));


%% Conservative estimate

if(CONSERVATIVE)    
    l = @(xp)( min( sqrt(sum((proj_X(xp) - xp).^2)),1) );
    E = @(x,xp)(  v(x) - l(xp) - v(proj_X(xp))  );
    E_num = E(x_samp_valid,x_samp_valid_plus);
    Emax = max(E_num);
    v = @(x)(cv'*basisEval(x) - Emax / (1-alpha));
end



%% Plot

% Load true MPI set
load('MPI_Julia_True','ZZ','grid')

% Create grid
[X1,X2] = meshgrid(linspace(-1,1,grid),linspace(-1,1,grid));
X12 = [X2(:)' ; X1(:)'];

% Evaluate -v + 1 on a grid (unit superlevel set then gives outer approximation)
Z = max(-v(X12)+1,0);
Z = reshape(Z,size(X1))';
Z(X1.^2 + X2.^2 > 1) = NaN;

% Plot
[~,h] = contourf(X1,X2,(Z>=1) + 2*ZZ  -3*double(~((Z < 1) & (ZZ >= 1)))  ); hold on
set(h,'LineColor','none')
colormap([ [1 1 1]; 0.85*[1 1 1]; [1 1 1] ;0.3*[1 1 1]  ; [1 0.4 0];   ]);
contour(X1,X2,Z,[1 1],'color','black','linewidth',2);hold on
plot(cos([0:0.01:2*pi]),sin([0:0.01:2*pi]),'--k','linewidth',2)

 
rmpath(genpath('./Resources'));












