function runPendulumSwingup_FT 
% checkDependency('yalmip') 
checkDependency('spotless') 
% checkDependency('mosek')
clear
clc
close all 
 
try
    % load nominal Trajectory and feedback control
    load('data.mat')
%     error
catch
    methods = 'Taylor'; % General Taylor
    % gen nominal trajectory by using GPOPS-II
%     sys = Pendulumdyn();
%     A = sys.dynamicA(0, [pi;0], 0);
%     B = sys.dynamicB(0, [pi;0], 0);
%     Q = diag([1 1]); R = 1;
%     [K,ST] = lqr(full(A),full(B),Q,R); 
    load('traj.mat') 
    figure
    hold on
    h = plot(xtraj);
    set(h,'LineWidth',2)
%     t = msspoly('t');
    x = msspoly('x',2);
    u = msspoly('u' ); 
%     xdottraj = fnder(xtraj);
    Q = diag([1 1]); R = 1; Qf = 2*diag([1 1]);  %Qf =ST;
    Nstep = 50;
    sys = sys.timecalculation(Nstep);
    [K,~,~] = Gaotvlqr(sys,xtraj,utraj,Q,R,Qf);
    
    Jetsys = DAODE(x,u,5,0.5);
    switch methods
        case 'Taylor'
            Jetsys = Jetsys.TaylorMap(sys,xtraj,utraj,K);
        case 'General'
            Jetsys = Jetsys.GeneralMap(sys,xtraj,utraj,K,R);
    end
    tic
    Xploy = Jetsys.integrate(xtraj.tspan,Nstep); 
    toc
end

figure
hold on
h = plot(xtraj);
set(h,'LineWidth',2)
ts = linspace(xtraj.tspan(1),xtraj.tspan(2),Nstep);

%% fast shootting
R0 = 0.5;  % initial set radius 

%% Lebesgue Moments
V_Lebesgue = LebesgueMeas(ts,x,R0,Xploy,1);
%% plot funnel
options = struct();
options.plotdims = [1;2];
options.x0 = [0;0];
[h,B] = plotFunnel(x,V_Lebesgue,xtraj,ts,options);

%% gen real funnel
methods = 'MonteCarlo'; % MonteCarlo Fast convexhull
% figure
% hold on
switch methods
    case 'Fast'
        xini = Ini_Convex(R0);
        [shootingX,A] = FastShoot(ts,x,xtraj,Xploy,xini);
    case 'MonteCarlo'  
        fdsys = @(t,x)(sys.dynamics(t,x,utraj.eval(t)+K{1}.eval(t)*(x-xtraj.eval(t))+K{2}.eval(t)));  
        xini = Ini_Uniform(R0);
        [shootingX,A]= MonteCarloShoot(fdsys,ts,xini);
end 

pre = abs(A-B)/A;
disp(['Funnel mismatch:',num2str(pre*100),'%'])
end
function xini = Ini_Convex(R)
Num = 100;
xini = [];  
theta = linspace(0,2*pi,Num); 
xini = [R*cos(theta);R*sin(theta)]; 
end
function xini = Ini_Uniform(R)
Num=500;
xini = [];  
kk=1;
Randomp = haltonset(2);
Randomp = scramble(Randomp,'RR2');
X0 = net(Randomp,1.5*Num)';
while length(xini) < Num 
    kk = kk+1;
    x0 =  0.5-X0(:,kk);
    if x0(1)^2+x0(2)^2 > R^2
        continue
    end 
    xini = [xini,x0];
end 
% hold on
% plot(xini(1,:),xini(2,:),'.k')
end

function [shootingX,A] = FastShoot(ts,x,xtraj,Xploy,xini)
 
hold on
Nstep = length(ts);
shootingX = [];
kk = 1;
for kk = 1:length(xini)
    x0 = xini(:,kk); 
    sX = []; 
    sX = [sX, double(subs(Xploy(:,1),x,x0))];
    plotX = xtraj.eval(ts(1)) + sX;
    for j = 2:Nstep
        sX = [sX, double(subs(Xploy(:,j),x,x0))];
        X = xtraj.eval(ts(j)) + sX(:,end);
        plotX = [plotX,X];
    end 
    plot(plotX(1,:),plotX(2,:),'b')  
    shootingX = [shootingX;sX];
end
end

function [shootingX,A] = MonteCarloShoot(dyn,ts,xini) 
% hold on
shootingX = []; 
for kk = 1:length(xini)
    x0 = xini(:,kk);
    sol = ode45(dyn,[ts(1),ts(end)],x0);
    output = ODE45SolTrajectory(sol,[length(x0),length(ts)]);
    X = output.eval(ts);
%     plot(X(1,:),X(2,:),'k')
    shootingX = [shootingX;X];
end
A = 0;
% figure
% hold on
for i = length(ts)-1:-1:1
    X1 = shootingX(1:2:end,i); X2 = shootingX(1:2:end,i+1);
    Y1 = shootingX(2:2:end,i); Y2 = shootingX(2:2:end,i+1);
    X = [X1;X2];    Y = [Y1;Y2];
    k = convhull(X,Y);
    k1 = convhull(X1,Y1); 
%     fill3(X(k),Y(k),repmat(0,1,length(k)),.7*[1 1 1],'LineStyle','none')
%     plot3(X(k),Y(k),repmat(-.1,1,length(k)),'k','LineWidth',5)
%     plot3(X1(k1),Y1(k1),repmat(.1,1,length(k1)),'k','LineWidth',.5)
    A = A + polyarea(X(k),Y(k)) - polyarea(X1(k1),Y1(k1));  
end
A = polyarea(X1(k1),Y1(k1))+A;
end


function vt = LebesgueMeas(ts,x,R,Xploy,deg)

basis = monomials(x,0:deg);
M0 = zeros(length(basis));
M0(1:3,1:3) = diag([1,1/R^2,1/R^2]);
MOM{1} = M0;

Mom = pinv(MOM{1});
vt{1} = basis'*MOM{1}*basis;
Q = basis*basis';
[~,p0,M0,~] = decomp(Q);
l0 = M0\Mom(:);
% l = LebesgueMom([0,2], R, 2);

for i = 2:length(ts)
    V = cleanPow(subs(Q,x,Xploy(:,i)),4);
    [x,p,M,sz] = decomp(V);
    [~,Locb] = ismember(p0,p,'rows');
    l = zeros(length(p),1);
    l(Locb) = l0;
    Mom = reshape(M*l,sz);
    MOM{i} = pinv(Mom)+1e-4*eye(length(Mom));
    vt{i} = basis'*MOM{i}*basis;
end

end

function [h,B] = plotFunnel(x,V,xtraj,ts,options)
h=[];
B = 0;
for i = length(V)-1:-1:1
    Xadd = xtraj.eval(ts(i));  Xadd2 = xtraj.eval(ts(i+1));
    y = getProjection(x,V{i}-subs(V{i},x,[0;0]),options);
    xfun0 = y+repmat(Xadd,1,length(y));
    y = getProjection(x,V{i+1}-subs(V{i+1},x,[0;0]),options);
    xfun1 = y+repmat(Xadd2,1,length(y));
    
    xx = [xfun0,xfun1(:,end:-1:1)];
    k = convhull(xx(1,:),xx(2,:));
    h=[h;fill3(xx(1,k),xx(2,k),repmat(0,1,length(k)),.7*[1 1 1],'LineStyle','none')];
    % plot convhulls
    h=[h;plot3(xx(1,k),xx(2,k),repmat(-.1,1,length(k)),'Color',[0 0 0],'LineWidth',5)];
    h=[h;plot3(xfun0(1,:),xfun0(2,:),repmat(.1,1,size(xfun0,2)),'Color',[.5 .5 .5])];
    k1 = convhull(xfun0(1,:),xfun0(2,:));
    B = B + polyarea(xx(1,k),xx(2,k)) - polyarea(xfun0(1,k1),xfun0(2,k1));   
end
B = polyarea(xfun0(1,k1),xfun0(2,k1))+B;

% h=[h;plot3(xfun0(1,:),xfun0(2,:),repmat(.1,1,size(xfun0,2)),'k','LineWidth',2)];
xlabel('\theta(rad)')
ylabel('$\dot{\theta}(rad/s)$','Interpreter','latex')
set(gcf,'color','w')
set(gca,'FontSize',18)
end