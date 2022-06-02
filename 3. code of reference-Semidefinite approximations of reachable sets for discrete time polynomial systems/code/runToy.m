function runToy
checkDependency('yalmip') 
checkDependency('mosek')
% Toy example by SDP-rexlation approximation
clear 
close all
clc
 
nX = 2; d = 8*2; R=1;T=100;e=0;
f = @(x)[0.5*(x(1)+2*x(1).*x(2));
     0.5*(x(2)-2*x(1).^3) ]; 
% Variables
xs = sdpvar(nX,1); 
fx = f(xs);dw = d+degree(fx);
[w,cw,wd] = polynomial(xs,d+degree(fx));
[v,cv] = polynomial(xs,d);
 
gX = 1 - xs'*xs; 
gX0 = 0.25^2-(xs(1)-0.5).^2-(xs(2)-0.5).^2; 
pow = degree(wd,xs,1);
l = LebesgueMom([0,d+degree(fx)], R, nX,pow);
 
u = sdpvar(1); 
[con1,L1,c1] = addSOScon_Yalmip(xs,v,gX0,d-2);
[con2,L2,c2] = addSOScon_Yalmip(xs,w-1-v,gX,d-2);
[con3,L3,c3] = addSOScon_Yalmip(xs,w,gX,d-2);
v1 = replace(v,xs,fx);
[con4,L4,c4] = addSOScon_Yalmip(xs,u+v1-v,gX,d-2); 
% con5 = [u>=0];

% Objective
obj = (cw'*l+T*u*l(1));  
Con = [con1;con2;con3;con4;sos(L1);sos(L2);sos(L3);sos(L4)]; 
coefs = [cw;cv;c1(:);c2(:);c3(:);c4(:);u]; 
% Solver parameteres
SDPsolver = 'mosek'; % or sedumi mosek
% options = getSolverParams(SDPsolver);
% Solve
tic
diagnostics = solvesos(Con,obj,sdpsettings('solver','mosek'),coefs);
 toc
Xs0= 1-2*rand(2,10000);
ind = find(0.25^2-(Xs0(1,:)-0.5).^2-(Xs0(2,:)-0.5).^2>=0);
Xs = Xs0(:,ind);
figure
hold on 
X = sdpvar(1,1); Y = sdpvar(1,1);
wv = monolist([X;Y],dw);
p = vectorize(sdisplay(double(cw)'*wv));
[X,Y] = meshgrid(linspace(-1.1,1.1,100),linspace(-1.1,1.1,100));
W = eval(p);
varw = max(max(W))-min(min(W));
disp(['variation = ' num2str(varw)]);
if varw > 1e-8
[c,h] = contourf(X,Y,e + W-1,[0 0],'color','k','linewidth',3);
hold on;
caxis([-0.2 1]);
colormap gray;
cm=colormap;
colormap(flipud(cm));
% [c,h] = contour(X,Y,e + W-1,[0 0],'color','k','linewidth',3);
end

plotCircle([0;0],1);
plot(Xs(1,:),Xs(2,:),'.','MarkerSize',10); 
for i=1:T
    Xs = dyn(Xs); 
    plot(Xs(1,:),Xs(2,:),'.','MarkerSize',10); 
end

title(['2r=',num2str(d/2)])
xlabel('x_1')
ylabel('x_2')
end

function X = dyn(x) 
X = [0.5*(x(1,:)+2*x(1,:).*x(2,:));
     0.5*(x(2,:)-2*x(1,:).^3) ]; 
end


