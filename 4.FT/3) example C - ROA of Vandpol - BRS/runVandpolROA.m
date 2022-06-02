function runVandpolROA
checkDependency('yalmip') 
checkDependency('spotless') 
checkDependency('mosek')

clear
clc
close all
dbstop if error
figure
hold on
axis([-5 5 -5 5])
[x1,x2] = meshgrid(-5:0.5:5,-5:0.5:5);
u = -x2;
v = -x2.*(1-x1.^2)+x1;
sumuv = sqrt(u.^2+v.^2);
u = u./sumuv;
v = v./sumuv;
quiver(x1,x2,u,v,0.5)
xlabel(gca,'x_1');
ylabel(gca,'x_2');
set(gca,'FontSize',18)
%% =========================real ROA=========================
sys = vandpol(); 
EndTime = 15;
vandpolw = @(t,x)(sys.dynamics(t,x,0));
% backwards
sol = ode45(vandpolw,[EndTime,0],[-2.0148;-0.0541]);
ts = EndTime:-0.1:0;
output = ODE45SolTrajectory(sol,[2,length(ts)]);
drawx = output.eval(ts);
y0 = drawx;
y0(:,68:end) = []; y0(:,68) = y0(:,1); 
h0 = plot(y0(1,:),y0(2,:),'k','lineWidth',2); 
%% ===================1.occupation measure ROA=====================
[y1,h1] = vandpolROA;
A = polyarea(y0(1,:),y0(2,:));
B = polyarea(y1(1,:)',y1(2,:)');
pre1 = abs(A-B)/A;
%% =====================2.bilienar search ROA=====================
t = msspoly('t');
u = msspoly('u',1);
x = msspoly('x',2); 
S = [0.378 -0.137;-0.137 0.278]*3;
V0 = V_function(t,x,S,[]);
V0.x0 = [0;0];
options.max_iterations = 20;
options.degV = 8;
options.converged_tol = 1e-2; 
options.method = {'bilinear'};
tic
[V,rho] = regionOfAttraction(sys,V0,options);
toc
y2 = getLevelSet(x,V,options);
h2 = plot3(y2(1,:),y2(2,:),repmat(3,1,size(y2,2)),'g','LineStyle','-','LineWidth',2);
 
B = polyarea(y2(1,:)',y2(2,:)');
pre2 = abs(A-B)/A;
%% =====================3. level set ROA=====================
options.degV = 8;
options.converged_tol = 1e-2; 
options.method = {'levelset'};
tic
[V,rho] = regionOfAttraction(sys,V0,options);
toc
y3 = getLevelSet(x,V/rho,options);
h3 =plot3(y3(1,:),y3(2,:),repmat(3,1,size(y3,2)),'c','LineStyle','-','LineWidth',2); 
B = polyarea(y3(1,:)',y3(2,:)');
pre3 = abs(A-B)/A;
%% FT ROA
try
    load('ex3_y.mat')
catch
    runVandpolROA_FT;
    load('ex3_y.mat')
end
h4 = plot(y(1,:),y(2,:),'r','lineWidth',2);
legend([h0,h1,h2,h3,h4],'real','occupation measure','bilinear',...
    'level-set','FT') 
B = polyarea(y(1,:)',y(2,:)');
pre4 = abs(A-B)/A;
 
disp(['ROA mismatch by occupation measure:',num2str(pre1*100),'%'])
disp(['ROA mismatch by bilienar search:',num2str(pre2*100),'%'])
disp(['ROA mismatch by level set:',num2str(pre3*100),'%'])
disp(['ROA mismatch by FT:',num2str(pre4*100),'%'])
end

function [Data1,HandleROA] = vandpolROA 
% compute ROA via measure
SDPsolver = 'mosek'; % or sedumi mosek 
d = 8;
tic
% Variables
x = sdpvar(2,1);
t = sdpvar(1,1);

% Define polynomials w(x) and v(t,x)
[w,cw] = polynomial(x,d);
[v,cv] = polynomial([t;x],d);
 
xb = 5;  
T = 100;   
f = [-x(2) ; -x(2)*(1-x(1)^2)+x(1) ] * T; 
  
gxT = (0.1^2 - x'*x);  
gX = [xb^2 - x(1)^2;xb^2 - x(2)^2];
gT = t*(1-t);
% Lebesgue moments on X
l = getLebesgueMoments(d,[-xb -xb ; xb xb],1);

Lv = jacobian(v,t) + jacobian(v,x)*f;
dk = d-2;
% Constraints (Note that the dynamics was scaled by T, so there T = 1)
[con1,L1,c1] = addSOScon_Yalmip([x;t],-Lv,[gX;gT],dk);
[con2,L2,c2] = addSOScon_Yalmip(x,replace(v,t,1),gxT,dk);
[con3,L3,c3] = addSOScon_Yalmip(x,w-replace(v,t,0)-1,gX,dk);
[con4,L4,c4] = addSOScon_Yalmip(x,w,gX,dk);
% Objective
obj = cw'*l; % minimization of int w d_lambda 
Con = [con1;con2;con3;con4;sos(L1);sos(L2);sos(L3);sos(L4)]; 
coefs = [cw;cv;c1(:);c2(:);c3(:);c4(:)]; 
% Solver parameteres
options = getSolverParams(SDPsolver);
if (~isempty(options))
    mset(options);
end

% Solve

diagnostics = solvesos(Con,obj,options,coefs);
toc
% Retrieve coefficients of w and V
cw = double(cw);
cv = double(cv);


%% Plots 
% Level-set plot of { x : v(0,x) >= 0 } and the true ROA
hold on
X = sdpvar(1,1); 
Y = sdpvar(1,1);
vv = monolist([t;X;Y],d);
p = vectorize(sdisplay(replace(cv'*vv,t,0) + 1)); % v(0,x) + 1
[X,Y] = meshgrid(-xb:0.05:xb,-xb:0.05:xb);
Z = eval(p);
[Data1, HandleROA] = contour(X,Y,Z, [1 1], '-b', 'linewidth',2);
 
end



