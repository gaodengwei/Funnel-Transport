
% Projects on the unit box
% Box = [lb; ub]
function xp = projBox(x,Box)

if(~exist('Box','var') || isempty(Box))
    Box = [-ones(1,size(x,1));ones(1,size(x,1))];
end

xp = zeros(size(x,1),size(x,2));
for i = 1:size(x,2)
    xp(:,i) = max(min(x(:,i),Box(2,:)'),Box(1,:)');
end


end