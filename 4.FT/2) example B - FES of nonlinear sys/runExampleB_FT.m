function runExampleB_FT
checkDependency('spotless') 
% checkDependency('mosek')
clear
clc
close all 

x = msspoly('x',1);
Nsample = 10000;
domain = [-1,1];
xs = linspace(domain(1),domain(2),Nsample);
%%  plot x^+
TransOrig = 2.1143+1.8*x+1.3571*x^2;
tic
Xnest = msubs(TransOrig,x,xs);
toc 

for d=1:10 
    GenKernel(x,d,xs,TransOrig,Xnest,domain,Nsample)
end

end

function ls = ProjectedGD(x,ker0,xs)

Y = msubs(ker0,x,xs);
ls = max(full(Y));
end

function GenKernel(x,d,xs,TransOrig,Xnest,domain,Nsample)

v = monomials(x,0:d);
Mx0 = v*v';
[~,p,C0,sz] = decomp(Mx0);
lv = LebesgueMomBox(domain',p);
M0 = reshape(C0*lv,sz);
L0 = chol(M0,'lower');
D0 = L0^-1; P0 = D0*v; invM0 = D0'*D0;
ker0 = P0'*P0;  % invM0 = M0^-1; ker0 = v'*invM0*v;

ls0 = ProjectedGD(x,ker0,xs);
%% trans
tic
Mx1 = subs(Mx0,x,TransOrig);
[~,p,C1,sz] = decomp(Mx1);
toc
lw = LebesgueMomBox(domain',p);
M1 = reshape(C1*lw,sz)+1e-10*eye(sz);
tic
invM1 = M1^-1;
ker1 = v'*invM1*v;
toc
ls1 = ProjectedGD(x,ker1,Xnest);
Display = [min(Xnest);max(Xnest)];DD = Display(2)-Display(1);
Xs = linspace(Display(1)-DD,Display(2)+DD/2,Nsample);

% plot kernel
figure
hold on
Y = msubs(ker1,x,Xs);
% Y0 = msubs(ker0,x,xs);
plot(Xs,Y,'k')
% plot(xs,Y0,'r')
xlabel('x')
ylabel('\rho')
switch d
    case 1
        title(['$\kappa_{t_0,1} \leq$' num2str(ls0) ',  ' '$\kappa_{t_2,1} \leq$' num2str(ls1)],'Interpreter','latex')
    case 2
        title(['$\kappa_{t_0,2} \leq$' num2str(ls0) ',  ' '$\kappa_{t_2,2} \leq$' num2str(ls1)],'Interpreter','latex')
    case 3
        title(['$\kappa_{t_0,3} \leq$' num2str(ls0) ',  ' '$\kappa_{t_2,3} \leq$' num2str(ls1)],'Interpreter','latex')
    case 4
        title(['$\kappa_{t_0,4} \leq$' num2str(ls0) ',  ' '$\kappa_{t_2,4} \leq$' num2str(ls1)],'Interpreter','latex')
    case 5
        title(['$\kappa_{t_0,5} \leq$' num2str(ls0) ',  ' '$\kappa_{t_2,5} \leq$' num2str(ls1)],'Interpreter','latex')
    case 6
        title(['$\kappa_{t_0,6} \leq$' num2str(ls0) ',  ' '$\kappa_{t_2,6} \leq$' num2str(ls1)],'Interpreter','latex')
    case 7
        title(['$\kappa_{t_0,7} \leq$' num2str(ls0) ',  ' '$\kappa_{t_2,7} \leq$' num2str(ls1)],'Interpreter','latex')
    case 8
        title(['$\kappa_{t_0,8} \leq$' num2str(ls0) ',  ' '$\kappa_{t_2,8} \leq$' num2str(ls1)],'Interpreter','latex')
    case 9
        title(['$\kappa_{t_0,9} \leq$' num2str(ls0) ',  ' '$\kappa_{t_2,9} \leq$' num2str(ls1)],'Interpreter','latex')
    case 10
        title(['$\kappa_{t_0,10} \leq$' num2str(ls0) ',  ' '$\kappa_{t_2,10} \leq$' num2str(ls1)],'Interpreter','latex')
end

h1 = plot(Xs,ls0*ones(1,Nsample),'r','LineWidth',2) ;
h2 = plot(Xs,ls1*ones(1,Nsample),'b--','LineWidth',2) ;
set(gca,'XLim',[-1 7]);
set(gca,'YLim',[0 ls1]);
Xmax = max(Xnest);
Xmin = min(Xnest);
Ymin = 0;
Ymax = ls1;
plot(Xmax*ones(1,2),[Ymin, Ymax],'k','LineWidth',2)
plot(Xmin*ones(1,2),[Ymin, Ymax],'k','LineWidth',2)
set(gca,'FontSize',18) 
fun = @(a)double(subs(ker1-18,x,a));
% x1 = fsolve(fun,1.5);
% x2 = fsolve(fun,5);
end

