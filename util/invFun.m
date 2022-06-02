function MOM = invFun(sx,p,xs)
Num = size(sx,1);     % sample num
Ps = matlabFunction(p(2:end)','Vars',xs);
Q = Ps(sx(:,1)',sx(:,2)');
Q = [ones(1,Num);Q];
d = 1/Num;
Mom = sum(Q,2)*d;
M = reshape(Mom,size(p));
[U,S,V] = svd(M);
Sdiag = diag(S);
Sdiag = max(Sdiag,1e-10);
M = U*(diag(Sdiag))*V';
MOM = inv(M);
end