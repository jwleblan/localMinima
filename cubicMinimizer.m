% cubicMinimizer - Finds the minima of a cubic polynomial
%{
%-------------------------------------------------------------------------------
% SYNTAX:
%   [result,a] = cubicMinimizer(x1,f1,g1,x2,f2,g2)
%
% DESCRIPTION: 
%   Finds the minimizer of the cubic through x1 and x2 via their function
%   values and derivatives (f & g).  These values may be vectors.
%
% INPUT:
%   x1,x2   - Distinct ordinate values
%   f1,f2   - Function value at x1,x2 respectively
%   g1,g2   - Function gradient at x1,x2 respectively
%
% RETURN VALUES:
%   result  - Local minima of the interpolant or NaN if no such minimia exists
%   a       - Coefficient on the cubic term
% 
% TODO:
%   Make other polynomial coefficients available.  Presently not neeeded
%
% Copyright (C) 2013 Joel W. LeBlanc
%-------------------------------------------------------------------------------
%}
function [result,a] = cubicMinimizer(x1,f1,g1,x2,f2,g2)
result = nan(numel(x1),1);
thresh = 1E-6;
dx = x2 - x1;
df = f2 - f1;
d1 = g1 + g2 - 3*df./dx;
d2 = max(0,d1.^2-g1.*g2);

if(nargout>1)
    a = (d1 + df./dx)./dx.^2;
end

ind = find(d2>=thresh);
d2 = sign(dx(ind)).*sqrt(d2(ind));
temp = (2*d2-g1(ind)+g2(ind));

ind2 = find(abs(temp) >= thresh);
ind = ind(ind2);

result(ind) = x2(ind) - dx(ind).*(d2(ind2)-d1(ind)+g2(ind))./temp(ind2);

end
