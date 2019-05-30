% quadraticMinimizer - Finds the minima of a quadratic polynomial
%{
%-------------------------------------------------------------------------------
% SYNTAX:
%   result = quadraticMinimizer(x1,f1,g1,x2,f2,g2)
%
% DESCRIPTION: 
%   Finds the minimizer of the quadratic through x1 and x2 via their function
%   values and derivatives (f & g).  Exactly one of f2 or g2 must be NaN.  If
%   no minimizer exists, NaN is returned.
%
% RETURN VALUES:
%   result  - Minimizer of the interpolant
%
% Copyright (C) 2013 Joel W. LeBlanc
%-------------------------------------------------------------------------------
%}
function result = quadraticMinimizer(x1,f1,g1,x2,f2,g2)
result = NaN(numel(x1),1);
dx = x2 - x1;
g1dx = g1.*dx;

ind1 = isnan(f2);
ind2 = find(~ind1);
ind1 = find(ind1);

% One function value and two derivatives
g1Ind = g1(ind1);
temp = (g2(ind1)-g1Ind)./dx(ind1);
ind3 = find(temp>0);
ind = ind1(ind3);
result(ind) = x1(ind)-g1Ind(ind3)./temp(ind3);

% Two function values and one derivative
temp = g1dx(ind2)+f1(ind2)-f2(ind2);
ind3 = find(temp<0);
ind = ind2(ind3);
result(ind) = x1(ind) + .5*g1dx(ind).*dx(ind)./temp(ind3);

end

% if(isnan(f2))       % One function value and two derivatives
%     temp = g2-g1;
%     if(temp==0)
% %         if(dx*g1 < 0)
% %             result = x2;
% %         else
% %             result = x1;
% %         end
%         result = NaN;
%     else
%         result = x1 - g1*dx/temp;
%     end
%     
% else                % Two function values and one derivative
%     temp = g1*dx+f1-f2;
%     if(temp==0)
% %         if(f1<f2)
% %             result = x1;
% %         else
% %             result = x2;
% %         end
%         result = NaN;
%     else
%         result = x1 + .5*g1*dx^2/temp;
%     end
%     
% end