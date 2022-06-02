function runVandpolROA_FT 
% checkDependency('yalmip') 
checkDependency('spotless') 
% checkDependency('mosek')
clear
clc
close all 

sys = vandpol(); 
EndTime = 15;

x = msspoly('x',2);
u = msspoly('u' );
Jetsys = DAODE(x,u,8,1.4,2);

Jetsys = Jetsys.GeneralMap(sys,[0;0],0);
try 
    % load('datavanplo.mat')   % 4-degree 100 steps
    % load('datavanplo2.mat')  % 4-degree 1500 steps
    % load('datavanplo3.mat')  % 8-degree 1500 steps
    load('datavanplo4.mat')    % 8-degree 1500 steps sample methods
catch
    Nstep = 20;
    Xploy = Jetsys.integrate([EndTime,0],Nstep);%flip(xtraj.tspan)
end 

vandpolw = @(t,x)(sys.dynamics(t,x,0));
sol = ode45(vandpolw,[EndTime,0],[-2.0148;-0.0541]);
ts = EndTime:-0.1:0;
output = ODE45SolTrajectory(sol,[2,length(ts)]);
drawx = output.eval(ts);
y0 = drawx;
y0(:,68:end) = []; y0(:,68) = y0(:,1);
%% ===============plot ROA===============
figure
hold on  
plot(y0(1,:),y0(2,:),'k','lineWidth',2); 
%% BRSs at different time slices 
R = 1.4;deg = 4;
color = {'r','y','g','cyan','k','b'};
hh=[];
kk=1;
for i=[1,100,200,500,1000,1500]
    [~,y] = LebesgueMeas(x,R,Xploy,deg,i);
    hold on
    hh{kk} = plot(y(1,:),y(2,:),color{kk},'lineWidth',1);
    kk=kk+1;
end
xlabel('x_1')
ylabel('x_2')   
set(gca,'FontSize',18)
axis([-2.5 2.5 -3 3])

legend([hh{1},hh{2},hh{3},hh{4},hh{5},hh{6}],'t=0s','t=-1s','t=-2s','t=-5s','t=-10s','t=-15s')
keyboard
save ex3_y y
end

function [V,y] = LebesgueMeas(x,R,Xploy,deg,i)
Levelset = 1.4;
basis = monomials(x,0:deg);
Q = basis*basis';
vt = subs(Q,x,clean(Xploy(:,i)));
[x,p,M,sz] = decomp(vt);
l = LebesgueMom(0, R, 2, p);
Mom = reshape(M*l,sz)/Levelset; 
V =  basis'*inv(Mom)*basis;  
options.x0 = [0;0];
y = getLevelSet(x,V-subs(V,x,[0;0]),options); 


end








