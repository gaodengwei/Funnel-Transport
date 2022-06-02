function moments = getLebesgueMoments(d, box, Basis, p)
% getLebesgueMoments Computes the moments up to degree d of the n-dimensional
% Lebesgue measure over a box. 
% INPUTS: 
%   "d" - degree. 
%   "box" - defined as [vector of lower bounds; vector of upper bounds] 
%   "Basis" - monomial basis  
%   "p" - (default=1) the p-norm to be used, should be 1 <= p < oo. 
% OUTPUTS: 
%   "moments" - vector of moments, as a column vector. 
%==========================================================================

if (~exist('p', 'var') || isempty(p))
    p = 1;
end
if (~exist('Basis', 'var') || isempty(Basis))
    Basis = 0;
end

n = size(box,2);

if (Basis == 1)
    disp('Generating moments in Yalmip basis')
    dv = monpowers(n,d);
elseif Basis==2
    disp('Generating moments in Gloptipoly basis')
    dv = genPowGlopti(n,d);
elseif Basis==3
    disp('Generating moments in msspoly basis')
    dv = mss_asd(n,0:d); 
%     monomials
else
    dv = Basis;
end

% there is one moment for each row of dv
moments = zeros(size(dv,1),1);

if (p==1)
    for i = 1:numel(moments)
        moments(i) = prod((box(2,:).^(dv(i,:)+1) - box(1,:).^(dv(i,:)+1)) ./ (dv(i,:)+1));
    end
    
elseif (1 < p) && (p < inf)
    if any(box < 0)
        warning('Negative domain main produce imaginary moments.')
    end
    
    for i = 1:numel(moments)
        moments(i) = prod((box(2,:).^(p*dv(i,:)+1) - box(1,:).^(p*dv(i,:)+1)) ./ (p*dv(i,:)+1));
        moments(i) = moments(i)^(1/p);
    end
     
elseif (p == inf)
    error('p == inf not implemented.');
    
else
    error('p not understood.');
end

end


