function [Con,L,c,v,Leq,veq] = addSOScon_Yalmip(x,p,g,d,varargin)
% @program extraction the "p" in semialgebraic constraints "g" with "x" in
% d-degree polynomial
% x --  Indeterminate variable
% p --  ploy in x
% g   --  semialgebraic constraints
% d   -- scalar integer, d > deg(g).

% [L,c,v] = polynomial(x,d-2);
if  nargin > 4
    geq = varargin{1}; 
end
if nargin > 5
    xeq = varargin{2};
else
    xeq  = x;
end
L = [];
C = [];
n = size(p,1);
m = size(g, 1); % total number of h's
v = monolist(x,d);

for i=1:n
    c = sdpvar(length(v),m);
    % Multipliers
    l = c'*v;
    L = [L,l];
    C = [C;c(:)];
end
if nargin < 5 && nargout<5
    Con = sos(p - L'*g); % ensure that p is SOS on g
else
    veq = monolist(xeq,d);
    Leq = [];
    meq = size(geq, 1); % total number of h's
    for i=1:n
        c = sdpvar(length(veq),meq);
        % Multipliers
        l = c'*veq;
        Leq = [Leq,l];
        C = [C;c(:)];
    end
    Con = sos(p - L'*g - Leq'*geq); % ensure that p is SOS on g
end
end


