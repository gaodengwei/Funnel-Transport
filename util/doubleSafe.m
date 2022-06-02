function y=doubleSafe(x)
y=double(x);
if (~isa(y,'double')) error('double failed'); end
end