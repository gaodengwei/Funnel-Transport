function runExampleA_SOS
checkDependency('spotless') 
checkDependency('mosek')
clear all
clc
close all 

x = msspoly('x',2);
A = [1 -0.5;0.1 2];
f = A*x; 
V = 0.01*x'*eye(2)*x;
C = expm(A*1);
xp = [[-1,1,1,-1,-1]' [-1, -1,1,1,-1]'];
xp = (1+C*xp')'; 
% gen  Semi-algebraic Set
G0 = [1-x;1+x];
options.solver = @spot_mosek; % spot_mosek spot_sedumi
options.degV = 2;
options.clean_tol = 1e-6;
options.converged_tol = .01;
options.degL1 = 6;
options.max_iterations = 10;
figure
hold on 
plot(xp(:,1),xp(:,2),'k','LineWidth',2)  
options.x0 = [0;0];
rho = 1; 
for iter=1
    [V,P,basis] = optimizeRho(x,G0,options);
    basis1 = subss(basis,x,inv(C)*x);
    V = clean(basis1'*P*basis1); 
    y = getLevelSet(x,subs(V,x,1-x),options);
    plot(y(1,:),y(2,:),'LineWidth',2)
end
A = polyarea(xp(:,1),xp(:,2));
B = polyarea(y(1,:)',y(2,:)');
pre = abs(A-B)/A;

disp(['FRS mismatch by SOS:',num2str(pre*100),'%'])
end
 

function [Vk,P,basis] = optimizeRho(x,G0,options)
disp('Step 2: Searching for rho.')
% main
prog = spotsosprog;
prog = prog.withIndeterminate(x);
Lmonom = monomials(x,0:options.degL1);
basis = monomials(x,0:options.degV/2);
[prog,P] = prog.newPSD(length(basis));
prog = prog.withPSD(P);
V = basis'*P*basis;
L0 = [];
for j=1:length(G0)
    [prog,L00] = prog.newFreePoly(Lmonom);
    L0 = [L0;L00];
end 
rho=1;
prog = prog.withSOS(L0);
prog = prog.withSOS(rho-V-L0'*G0);   
[prog,obj] = maxdet(prog,P);
% Solve SOS program
pars = spot_sdp_default_options();
sol = prog.minimize(-obj,options.solver,pars); 
P = double(sol.eval(P));
Vk = clean(sol.eval(V));
end



