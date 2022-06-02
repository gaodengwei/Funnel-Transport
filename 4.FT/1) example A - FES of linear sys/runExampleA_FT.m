function runExampleA_FT()
checkDependency('spotless') 
% checkDependency('mosek')
clear
clc
close all 

% inverse moment matrix in a square
xp = [[-1,1,1,-1,-1]' [-1, -1,1,1,-1]'];

x = msspoly('x',2);
basis = monomials(x,0:4);
P = basis*basis'; 
[~,p,~,~] = decomp(P); 
l = LebesgueMomBox([-1 -1; 1 1],p); 
%% solve moment matrix in a getLebesgueMoments 
levelSet = LebesgueMoments(x,basis,P,l,xp);
clf
hold on
%% affine trans in Lebesgue measure 
A = [1 -0.5;0.1 2];
C = expm(A*1);
% C = [2.6211 -2.3163; 0.4633 7.253];
y = 1+C*x;
xp = (1+C*xp')';
tic
P = subs(P,x,y);
toc
[~,~,M0,sz] = decomp(P); 
M = reshape(M0*l,sz); 
tic
V = basis'*inv(M)*basis;
toc
V0 = clean(V); 
V = V0-levelSet+1; 
options.x0 = [1;1];
y = getLevelSet(x,V,options);
 
gcf;
hold on   
plot(y(1,:),y(2,:),'r','LineWidth',2)   
% plot(1,1,'*')
plot(xp(:,1),xp(:,2),'LineWidth',2,'Color',[0,0,0])

A = polyarea(xp(:,1),xp(:,2));
B = polyarea(y(1,:)',y(2,:)');
pre = abs(A-B)/A;

disp(['FRS mismatch by FT:',num2str(pre*100),'%'])

end

function levelSet = LebesgueMoments(x,basis,P,l,xp)
[~,~,M0,sz] = decomp(P);
M = reshape(M0*l,sz); 
V = basis'*inv(M)*basis;
V0 = clean(V);
levelSet = double(subs(V0,x,[1;1]));
% V = V0-levelSet+1;
V = V0/levelSet;
figure 
hold on
xlabel('x_1')
ylabel('x_2')
y = getLevelSet(x,V); 
plot(y(1,:),y(2,:),'r','LineWidth',2) 
plot(0,0,'o')
plot(xp(:,1),xp(:,2),'LineWidth',2,'Color',[0,0,0]) 

 
 

end
  


 




