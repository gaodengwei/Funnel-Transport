function a = cleanPow(a,tol,LSflag,xs,R)

if nargin < 3, LSflag = 0; end  %default
if nargin < 2, tol = 4; end % default

switch LSflag
    case 0
        a = DirecTruncation(a,tol);
    case 1
        a = QPapprox(a,tol,R);
    case 2 
        a = surfacefit(a,tol,xs); 
end
end
 

function a = DirecTruncation(a,tol)
% Remove terms with large order Directly
[x,p,M]=decomp(a);
if length(tol)==1
    M(:,sum(p,2)>tol) = 0;
    a=recomp(x,p,M,size(a));
else
    M(:,sum(p,2)>tol(2)) = 0;
    M(:,sum(p,2)<tol(1)) = 0;
    a=recomp(x,p,M,size(a));
end
end

function a = QPapprox(a,tol,R)
% Remove terms with large order with fminunc 
n = length(a);
[x,p] = decomp(a);
if all(sum(p,2)<tol), return; end

basis = monomials(x,0:tol);
Q = basis*basis';
[~,dd,M,szq] = decomp(Q);
ld = LebesgueMom(2, R, n, dd);

Mom =  reshape(M*ld,szq);
x0=ones(size(basis));
options.Display = 'off';
Poly = [];
for i=1:n
    poly = a(i)*basis;
    [~,p,M0] = decomp(poly);
    lb = M0*LebesgueMom(2, R,n, p);
    fun = @(a) a'*Mom*a-2*a'*lb;
    xopt = fminunc(fun,x0,options);
    Poly = [Poly;xopt'*basis];
end
a =  Poly;

end


function a = surfacefit(a,tol,xs)
% Remove terms with large order with sample points LS approx
[x,p,~,sz] = decomp(a);
if all(sum(p,2)<tol)
    return
end
load('sampleX.mat')
ps = msspoly2sym(x,xs,a);
pfun1 = plfcnchk(char(ps(1)));
pfun2 = plfcnchk(char(ps(2)));
v1 = feval(pfun1,X(:),Y(:));
v2 = feval(pfun2,X(:),Y(:));

dd = mss_asd(2,0:tol);
aa = [X(:)';Y(:)'];
P = nummonomials(aa,dd)';
c1 = P\v1;
c2 = P\v2;
M = [c1';c2'];
a = recomp(x, dd, M, sz);

end


function a = SDPapprox2(a,tol)
% Remove terms with large order with LS-SOS
% replaced by SDPapprox since the SOS is slow
R=0.1;
n = length(a);
[x,p,M0] = decomp(a);
if all(sum(p,2)<tol), return; end

prog = spotsosprog;
prog = prog.withIndeterminate(x);
basis = monomials(x,0:tol);

Q = basis*basis';
[~,dd,M,szq] = decomp(Q);
ld = LebesgueMom(2, R, n, dd);
Mom =  reshape(M*ld,szq);
% solve SOS
solver = @spot_mosek;
pars = spot_sdp_default_options();
pars.verbose = 0;
Poly= [];
for i=1:n
    poly = a(i)*basis;
    [~,p,M0] = decomp(poly);
    lb = M0*LebesgueMom(2, R, n, p);
    [prog,P,coeff] = prog.newFreePoly(basis);
    obj = coeff'*Mom*coeff-2*coeff'*lb;
    sol = prog.minimize(obj,solver,pars);
    Poly = [Poly;sol.eval(P)];
end
a= Poly;
end

function b = nummonomials(a,p)
b = [];
for i=1:length(p)
    c = (a(1,:).^p(i,1)).*(a(2,:).^p(i,2));
    b=[b;c];
end
end