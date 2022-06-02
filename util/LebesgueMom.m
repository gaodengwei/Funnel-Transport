function [l, alphas, theta] = LebesgueMom(d, R, nX, alphas)
%  Computes the moments in a ball with Radius @R
% See "How to Integrate a Polynomial over a Sphere", G. B. Folland for
% details

if nargin < 4
    alphas = mss_asd(nX,d(1):d(2));
end

betas = 0.5*(alphas + 1);
Ra = (R.^(sum(alphas,2) + nX))./(sum(alphas,2) + nX);
IS = 2*prod(gamma(betas),2)./(gamma(sum(betas,2)));
l = Ra.*IS;
alphaszero = any((mod(alphas,2) ~= 0),2);
% alphaszero = any(alphaszero,2);
l(alphaszero) = 0;

if nargout>2
    % not right here; need to be check
    wp = (2*pi^((nX+1)/2))/gamma((nX+1)/2);  
    d = sum(alphas,2);
    a = wp*((d+1).*(d+2).*(d+3))./((d+nX+1).*(d+nX+2).*(2*d+nX+6))./l(1);
    theta = max(prod(alphas,2))/a;   % level set of SOS poly
end


end



