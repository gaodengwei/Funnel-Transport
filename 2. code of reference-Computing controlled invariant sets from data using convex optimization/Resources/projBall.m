% x: n x N
function xp = projBall(x,r)

if(~exist('r','var') || isempty(r))
    r = 1;
end

normX = sqrt(sum(x.^2));
xp = r * x ./ repmat(max(normX,r),2,1);

end