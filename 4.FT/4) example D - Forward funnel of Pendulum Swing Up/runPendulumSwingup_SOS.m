function runPendulumSwingup_SOS 
checkDependency('yalmip') 
% checkDependency('spotless') 
checkDependency('mosek')
clear
clc
close all 
  
load('traj.mat') 
%% compete Forward funnel
Nstep = 50;
sys = sys.timecalculation(Nstep);
Q = diag([1 1]);
R = 1; 
Qf = 2*diag([1 1]);
[K,V] = Gaotvlqr(sys,xtraj,utraj,Q,R,Qf,[]);
% Vv = V.getPoly(0);
% S1 = double(.5*subs(diff(diff(Vv,x)',x),x,0*x));
% S0 = 4*eye(2);
% AA = S0/S1;
% V = V*AA;
% figure
% options.x0 = [0;0];
% y2 = getLevelSet(x,V.getPoly(0),options);
% plot(y2(1,:),y2(2,:))

options = struct();
options.degL1 = 2;
options.K = K;
ts = sys.breaks;
options.rho0_tau = 3; 
options.X0 =  0.5^2-x'*x; 
tic 
V = FRS_Funnel_SOS(sys,V,ts,options);
toc
keyboard
 

end
 