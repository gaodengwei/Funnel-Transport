function runToy_FT
checkDependency('spotless') 
% checkDependency('mosek')
% Toy example by FT  Dirac measure approximation
clear
close all
clc
dbstop if error

% nX = 2; r = 6; d = 6; R=0.25;
f = @(x)[0.5*(x(1)+2*x(1).*x(2));
     0.5*(x(2)-2*x(1).^3) ]; 
x0 = [0.5;0.5];
xf = x0;
Xf = xf;
for i=1:7
    xf = f(xf); 
    Xf = [Xf,xf];
end
 
x = msspoly('x',2);
  
%% transport
%% plot
XXf = x;
figure
hold on
plotCircle([0;0],1);
Xs0= 1-2*rand(2,10000);
ind = find(0.25^2-(Xs0(1,:)-0.5).^2-(Xs0(2,:)-0.5).^2>=0);
Xs = Xs0(:,ind);   
plot(Xs(1,:),Xs(2,:),'.','MarkerSize',10);   
[V_dirac,h2] = DiracMeas(x,Xs',4,Xf(:,1));
 
tic 
for i=1:7   
    Xs = dyn(Xs); 
    plot(Xs(1,:),Xs(2,:),'.','MarkerSize',10);   
    % Dirac measure 
     
    [V_dirac,h2] = DiracMeas(x,Xs',4,Xf(:,i+1));
end
toc
plot(Xf(1,:),Xf(2,:),'k','LineWidth',2)
xlabel('x_1')
ylabel('x_2')



 



end

function X = dyn(x)

X = [0.5*(x(1,:)+2*x(1,:).*x(2,:));
     0.5*(x(2,:)-2*x(1,:).^3) ];

end

function [basis,Mom,pow] = genMom(x,d,R,nX)
    basis = new_monomials(x,0:d);
    Q = basis*basis';
    [~,pow,coe,sz] = decomp(Q);
    l = LebesgueMom([0,d], R, nX, pow);
    Mom =  reshape(coe*l,sz); 
end

function [V,h2] = DiracMeas(x,xT,deg,x0)

basis = monomials(x,0:deg);
P = basis*basis';
xs = sym('x',[1 2],'real')';
ps = msspoly2sym(x,xs,P);
 

inverseM = invFun(xT,ps,xs); 
levelSet = 60;% 60 30  
hold on 
V = basis'*inverseM*basis; 
options.x0 = x0;
y = getLevelSet(x,V/levelSet,options);
h2 = plot(y(1,:),y(2,:),'r','LineWidth',2);
 
  
end