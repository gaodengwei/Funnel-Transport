% Evaluates the monomial basis (Yalmip ordering) on a given vector of points

% INPUT:
% d : degree
% X : n x N matrix of points on which to evaluate


% OUTPUT:
% monVal: nchoosek(n+d,d) x N matrix

function monVal = evalMonBasis(d,X)

n = size(X,1);
N = size(X,2);

mpow = monpowers(n,d);

monVal = zeros(size(mpow,1),N);

for i = 1:size(mpow,1)
    monVal(i,:) = prod(X.^repmat(mpow(i,:)',1,N),1);
end


end