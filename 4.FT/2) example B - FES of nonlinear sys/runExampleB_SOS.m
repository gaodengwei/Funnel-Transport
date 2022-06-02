function runExampleB_SOS
checkDependency('yalmip') 
checkDependency('mosek')
clear
clc
close all 

x = sdpvar(1,1);
y = sdpvar(1,1);
N = 1e4;
m = 1;
dd = 5; 
f = 2.125 + 1.5*x + 1.25*x^2 + 0.5*x^3 + 0.125*x^4;
lo = -10; up = 10;
g = [1 - x^2];

tic
[y,J] = exists2(x, y, f, dd, dd, g, lo, up, eps,m);
toc
hold on   
plot(y,0.2-J,'k')
plot(y,0.2*ones(1,100),'r','LineWidth',2)

xs = 1-2*rand(1,N); 
fx = 2.125 + 1.5*xs + 1.25*xs.^2 + 0.5*xs.^3 + 0.125*xs.^4;
plot(xs,fx,'.k','markersize',5);
axis([0,7,0.18,0.22])
plot([1.5,1.5],[0.18 0.22],'k','LineWidth',2)
plot([5,5],[0.18 0.22],'k','LineWidth',2)
xlabel('x')
ylabel('\rho')
set(gca,'FontSize',18)
end
function [y1, G] = exists2(x, yscale, f, k, order, g, lo, up, eps, n)
N = 100;
m = length(g);
order = max(k, order);
mons = monpowers(n,2 * k); 
[J, ck] = polynomial(yscale, 2 * k); 
gammak = getLebesgueMoments(2*k,[-1; 1],1); 
obj = ck' * gammak;
cstr = [];
coeffsos = [];
qk = 0;
for j = 1:m
    gj = g(j);
    dj = 2 * ceil(degree(gj)/2);
    [s,c] = polynomial([x; yscale], 2 * order - dj);
    coeffsos = [coeffsos; c];
    cstr = [cstr sos(s)];
    qk = qk + s * gj;
end
a = zeros(n, 1);
p = zeros(n, 1);
for i = 1:length(f)
    a(i) = 2  / (up(i) - lo(i));
    p(i) = - (up(i) + lo(i))/(up(i) - lo(i)); 
    fscale(i) = a(i) * f(i) + p(i); 
end
Ff = 0;
for i=1:n
    Ff = Ff + (fscale(i) - yscale(i))^2;
end
F = [sos(Ff - eps  - J - qk), cstr];
solvesos(F, -obj, sdpsettings('solver','mosek','sos.scale','1'), [coeffsos;ck]);
Jk = double(ck);

yscale1=ndgrid(linspace(-1,1,N)); 
y1 = (yscale1 - p(1))/a(1); 
G = 0;
for alpha = 1:length(mons)
    beta = mons(alpha,:);
    G = G + Jk(alpha).* yscale1.^beta;
end
for i = 1:1e2
    if 1 <= yscale1(i)^2
        G(i) = 2;
    end
end
end
 

