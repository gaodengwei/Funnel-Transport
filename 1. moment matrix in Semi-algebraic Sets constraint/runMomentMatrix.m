clc
clear all
close all
checkDependency('spotless') 
checkDependency('mosek')
x = msspoly('x',2);
%  Semi-algebraic Sets constraint
g = [1-x;x+1]; 
% basis and lagrangian multiplier 
basis = monomials(x,0:4);  
Lxmonom = monomials(x,0:4);
%%  =================inner approximation=================
prog = spotsosprog;
prog = prog.withIndeterminate(x);
[prog,M] = prog.newPSD(length(basis));  
[prog,L] = prog.newFreePoly(Lxmonom,length(g));
[prog,slack] = prog.newPos(1);

for i=1:4
    prog = prog.withSOS(L(i)*g(i)-(slack-basis'*M*basis));
end
prog = prog.withSOS(L);

options = spot_sdp_default_options(); 
prog = prog.withPos(.1-trace(M));
sol = prog.minimize(-slack,@spot_mosek,options); 
rho = double(sol.eval(slack)); 
V = basis'*sol.eval(M)*basis/rho;
options.x0 = [0;0];
y = getLevelSet(x,V,options);
figure
plot(y(1,:),y(2,:),'k','LineWidth',2)
hold on
plot([-1,1,1,-1,-1],[-1,-1,1,1,-1],'k--','LineWidth',2)
axis([-1.5 1.5 -1.5,1.5])
xlabel('x')
ylabel('y')

%% ================= outer approximation================= 
clear prog M L slack
prog = spotsosprog; 
prog = prog.withIndeterminate(x);
[prog,M] = prog.newPSD(length(basis));  
[prog,L] = prog.newFreePoly(Lxmonom,length(g));
[prog,slack] = prog.newPos(1);
for i=1:4
    prog = prog.withSOS(slack-basis'*M*basis-L(i)*g(i));
end
prog = prog.withSOS(L);

options = spot_sdp_default_options();
prog = prog.withPos(trace(M)-1.3e-6); 
sol = prog.minimize(slack,@spot_mosek,options);

rho = double(sol.eval(slack));
% rho =1;
V = basis'*sol.eval(M)*basis/rho;
options.x0 = [0;0];
y = getLevelSet(x,V,options);
figure
plot(y(1,:),y(2,:),'k','LineWidth',2)
hold on
plot([-1,1,1,-1,-1],[-1,-1,1,1,-1],'k--','LineWidth',2)
axis([-1.5 1.5 -1.5,1.5])
xlabel('x')
ylabel('y')
%% =================approximation=================
clear prog M L slack
prog = spotsosprog; 
prog = prog.withIndeterminate(x);
[prog,M] = prog.newPSD(length(basis));  
[prog,L] = prog.newFreePoly(Lxmonom,length(g));
[prog,slack] = prog.newPos(1);
for i=1:4
    prog = prog.withSOS(1-basis'*M*basis-L(i)*g(i));
    prog = prog.withSOS(L(i)*g(i)-(1-basis'*M*basis));
end
prog = prog.withSOS(L); 
options = spot_sdp_default_options(); 
sol = prog.minimize(0,@spot_mosek,options); 
rho =1;
V = basis'*sol.eval(M)*basis/rho;
options.x0 = [0;0];
y = getLevelSet(x,V,options);
figure
plot(y(1,:),y(2,:),'k','LineWidth',2)
hold on
plot([-1,1,1,-1,-1],[-1,-1,1,1,-1],'k--','LineWidth',2)
axis([-1.5 1.5 -1.5,1.5])
xlabel('x')
ylabel('y') 
