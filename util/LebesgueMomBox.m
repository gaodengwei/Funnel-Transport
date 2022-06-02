function moments = LebesgueMomBox(box, dv )
%  Computes the moments in a box 
    
    moments = zeros(size(dv,1),1);
    for i = 1:numel(moments)
        moments(i) = prod((box(2,:).^(dv(i,:)+1) - box(1,:).^(dv(i,:)+1)) ./ (dv(i,:)+1));
    end
end